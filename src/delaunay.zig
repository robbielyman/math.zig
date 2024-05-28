const std = @import("std");
const math = @import("math.zig");

const DelaunayMap = std.AutoArrayHashMap([2]usize, TriangleData);

pub const TriangleData = struct {
    vtx: ?usize,
    flippable: bool = true,
};

const ghost: usize = std.math.maxInt(usize);

/// supports constructing constrained Delaunay triangulations as explained by Shewchuk:
/// https://people.eecs.berkeley.edu/~jrs/meshpapers/delnotes.pdf
pub fn PlanarTriangulation(comptime F: type) type {
    return struct {
        const Self = @This();
        map: std.AutoArrayHashMap([2]usize, TriangleData),
        verts: std.ArrayList([2]F),
        allocator: std.mem.Allocator,

        /// returns a struct to which you can repeatedly call `addVertex` or `addEdge`
        /// to iteratively build up the triangulation.
        pub fn init(allocator: std.mem.Allocator) Self {
            return .{
                .map = std.AutoArrayHashMap([2]usize, TriangleData).init(allocator),
                .verts = std.ArrayList([2]F).init(allocator),
                .allocator = allocator,
            };
        }

        /// frees the memory occupied by the triangulation
        pub fn deinit(self: *Self) void {
            self.map.deinit();
            self.verts.deinit();
            self.* = undefined;
        }

        /// returns a list of indices into `self.verts.items`
        /// the length of the list will be a multiple of three unless `self.verts.items` is too short
        /// and each group of three within the list represents a constrained Delaunay triangle
        pub fn getNonGhostTriangles(self: *const Self, allocator: std.mem.Allocator) std.mem.Allocator.Error![]usize {
            var triangles = std.ArrayList(usize).init(allocator);
            defer triangles.deinit();
            if (self.verts.items.len < 3) {
                for (0..self.verts.items.len) |i| {
                    try triangles.append(i);
                }
                return try triangles.toOwnedSlice();
            }
            for (self.map.keys()) |k| {
                if (k[0] == ghost or k[1] == ghost) continue;
                const v = self.map.get(k).?.vtx orelse continue;
                if (v == ghost) continue;
                var idx: usize = 0;
                const t = .{ k[0], k[1], v };
                while (idx + 3 <= triangles.items.len) : (idx += 3) {
                    if (roteq(t, triangles.items[idx..][0..3].*)) break;
                } else try triangles.appendSlice(&t);
            }
            return try triangles.toOwnedSlice();
        }

        fn roteq(t: [3]usize, s: [3]usize) bool {
            const rots: [3][3]usize = .{
                .{ t[0], t[1], t[2] },
                .{ t[1], t[2], t[0] },
                .{ t[2], t[0], t[1] },
            };
            for (rots) |r| {
                if (std.mem.eql(usize, &r, &s)) return true;
            }
            return false;
        }

        /// adds an edge to the triangulation by first adding the vertices if they are not already included
        /// removing any edges that cross the constraining edge
        /// and retriangulating the interior using a "gift-wrapping" algorithm
        pub fn addEdge(self: *Self, edge: [2][2]F) (error{Clobber} || std.mem.Allocator.Error)!void {
            const e = try self.searchOrAddVerts(edge);
            var add = false;
            if (self.map.getPtr(.{ e[1], e[0] })) |ptr| {
                ptr.flippable = false;
                self.map.getPtr(e).?.flippable = false;
                add = ptr.vtx == null;
            } else add = true;
            if (!add) return;
            const upper_list, const lower_list = try self.clearIntersection(e);
            defer self.allocator.free(upper_list);
            defer self.allocator.free(lower_list);
            try self.giftWrap(upper_list, .forwards);
            try self.giftWrap(lower_list, .backwards);
        }

        fn searchOrAddVerts(self: *Self, edge: [2][2]F) (error{Clobber} || std.mem.Allocator.Error)![2]usize {
            var e: [2]usize = undefined;
            for (edge, 0..) |vtx, i| {
                    for (0..self.verts.items.len) |j| {
                        const idx = self.verts.items.len - 1 - j;
                        const v = self.verts.items[idx];
                        if (math.eq(v[0], vtx[0]) and math.eq(v[1], vtx[1])) {
                            e[i] = idx;
                            break;
                        }
                    } else {
                        try self.addVertex(vtx);
                        e[i] = self.verts.items.len - 1;
                    }
            }
            return e;
        }

        fn giftWrap(self: *Self, list: []const usize, comptime dir: enum { forwards, backwards }) (error{Clobber} || std.mem.Allocator.Error)!void {
            if (list.len < 3) return;
            switch (comptime dir) {
                .forwards => {
                    var idx: usize = 1;
                    for (2..list.len - 1) |i| {
                        if (self.inCircumcircle(.{ list[0], list[list.len - 1], list[idx] }, list[i]))
                            idx = i;
                    }
                    try self.addTriangle(.{ list[0], list[list.len - 1], list[idx] });

                    try self.giftWrap(list[0 .. idx + 1], .forwards);
                    try self.giftWrap(list[idx..list.len], .forwards);
                },
                .backwards => {
                    var idx: usize = 1;
                    for (2..list.len - 1) |i| {
                        if (self.inCircumcircle(.{ list[list.len - 1], list[0], list[idx] }, list[i]))
                            idx = i;
                    }
                    try self.addTriangle(.{ list[list.len - 1], list[0], list[idx] });
                    try self.giftWrap(list[0 .. idx + 1], .backwards);
                    try self.giftWrap(list[idx..list.len], .backwards);
                },
            }
        }

        fn clearIntersection(self: *Self, edge: [2]usize) (error{Clobber} || std.mem.Allocator.Error)!struct { []usize, []usize } {
            var upper_list = std.ArrayList(usize).init(self.allocator);
            defer upper_list.deinit();
            var lower_list = std.ArrayList(usize).init(self.allocator);
            defer lower_list.deinit();
            try upper_list.append(edge[0]);
            try lower_list.append(edge[0]);
            var seed, var v = seed: {
                for (self.map.keys()) |e| {
                    const v = self.map.get(e).?.vtx orelse continue;
                    if (v == ghost or e[1] == ghost) continue;
                    if (e[0] == edge[0]) break :seed .{ e, v };
                }
                unreachable;
            };
            while (!self.cuts(edge, .{ seed[1], v })) {
                seed[1] = v;
                v = self.map.get(seed).?.vtx.?;
            }
            var tri: [3]usize = .{ seed[0], seed[1], v };
            try upper_list.append(v);
            try lower_list.append(seed[1]);
            v = edge[0];
            while (tri[0] != edge[1] and tri[1] != edge[1] and tri[2] != edge[1]) {
                const edges: [3][2]usize = .{
                    .{ tri[2], tri[1] },
                    .{ tri[0], tri[2] },
                    .{ tri[1], tri[0] },
                };
                const tseq: [3]usize, const vseq = tseq: {
                    for (edges, 0..) |e, i| {
                        if (tri[i] == v) {
                            const vseq = self.map.get(e).?.vtx.?;
                            break :tseq .{ .{ e[0], e[1], vseq }, vseq };
                        }
                    }
                    unreachable;
                };
                if (self.inHalfspace(edge, vseq)) {
                    if (vseq != edge[1])
                        try upper_list.append(vseq);
                    v = if (self.inHalfspace(edge, tseq[0])) tseq[0] else tseq[1];
                } else {
                    if (vseq != edge[1])
                        try lower_list.append(vseq);
                    v = if (self.inHalfspace(edge, tseq[0])) tseq[1] else tseq[0];
                }
                self.deleteTriangle(tri) catch |err| switch (err) {
                    error.Clobber => return error.Clobber,
                    error.NoTriangle => unreachable,
                };
                tri = tseq;
            }
            self.deleteTriangle(tri) catch |err| switch (err) {
                error.Clobber => return error.Clobber,
                error.NoTriangle => unreachable,
            };
            try upper_list.append(edge[1]);
            try lower_list.append(edge[1]);
            return .{ try upper_list.toOwnedSlice(), try lower_list.toOwnedSlice() };
        }

        fn cuts(self: *const Self, ab: [2]usize, cd: [2]usize) bool {
            if (cd[0] == ghost or cd[1] == ghost) return false;
            if (self.inHalfspace(ab, cd[0]) == self.inHalfspace(ab, cd[1])) return false;
            const mat, const size = math.geometry.affineIntersection(2, 2, 2, F, .{
                self.verts.items[ab[0]],
                self.verts.items[ab[1]],
            }, .{
                self.verts.items[cd[0]],
                self.verts.items[cd[1]],
            });
            for (mat[0..size]) |vec| {
                if (between(vec, self.verts.items[ab[0]], self.verts.items[ab[1]])) return true;
            }
            return false;
        }

        /// inserts a vertex into the triangulation by using a version of the Bowyer--Watson algorithm
        pub fn addVertex(self: *Self, vtx: [2]F) (error{Clobber} || std.mem.Allocator.Error)!void {
            for (self.verts.items) |v| {
                if (math.eq(v[0], vtx[0]) and math.eq(v[1], vtx[1])) return;
            }
            try self.verts.append(vtx);
            const keys = self.verts.items;
            if (keys.len < 2) return;
            if (keys.len == 2) {
                try self.addTriangle(.{ 0, 1, ghost });
                try self.addTriangle(.{ 1, 0, ghost });
                return;
            }
            if (keys.len == 3) {
                self.deleteTriangle(.{ 0, 1, ghost }) catch |err| switch (err) {
                    error.Clobber => return error.Clobber,
                    error.NoTriangle => unreachable,
                };
                self.deleteTriangle(.{ 1, 0, ghost }) catch |err| switch (err) {
                    error.Clobber => return error.Clobber,
                    error.NoTriangle => unreachable,
                };
                const tri: [3]usize = switch (math.geometry.orient(2, F, .{
                    self.verts.items[0],
                    self.verts.items[1],
                    self.verts.items[2],
                })) {
                    .gt, .eq => .{ 0, 1, 2 },
                    .lt => .{ 0, 2, 1 },
                };
                try self.addTriangle(tri);
                try self.addTriangle(.{ tri[2], tri[1], ghost });
                try self.addTriangle(.{ tri[0], tri[2], ghost });
                try self.addTriangle(.{ tri[1], tri[0], ghost });
                return;
            }
            const tri = self.findTriangle(keys.len - 1);
            self.deleteTriangle(tri) catch |err| switch (err) {
                error.Clobber => return error.Clobber,
                error.NoTriangle => unreachable,
            };
            const edges: [3][2]usize = .{
                .{ tri[0], tri[1] },
                .{ tri[1], tri[2] },
                .{ tri[2], tri[0] },
            };
            for (edges) |e| {
                try self.digCavity(keys.len - 1, e);
            }
        }

        fn addTriangle(self: *Self, t: [3]usize) (error{Clobber} || std.mem.Allocator.Error)!void {
            const edges: [3][2]usize = .{
                .{ t[1], t[2] },
                .{ t[2], t[0] },
                .{ t[0], t[1] },
            };
            for (edges) |e| {
                const res = try self.map.getOrPut(e);
                if (res.found_existing) {
                    if (res.value_ptr.vtx != null)                         return error.Clobber;
                } else {
                    res.value_ptr.* = .{ .vtx = null };
                }
            }
            for (edges, 0..) |e, i| {
                const ptr = self.map.getPtr(e).?;
                ptr.vtx = t[i];
            }
        }

        fn deleteTriangle(self: *Self, t: [3]usize) error{ Clobber, NoTriangle }!void {
            const edges: [3][2]usize = .{
                .{ t[1], t[2] },
                .{ t[2], t[0] },
                .{ t[0], t[1] },
            };
            for (edges, 0..) |e, i| {
                const get = self.map.get(e) orelse return error.NoTriangle;
                const w = get.vtx orelse return error.NoTriangle;
                if (w != t[i]) return error.Clobber;
            }
            for (edges) |e| {
                const ptr = self.map.getPtr(e).?;
                ptr.vtx = null;
            }
        }

        fn digCavity(self: *Self, v: usize, edge: [2]usize) (error{Clobber} || std.mem.Allocator.Error)!void {
            if (edge[0] == v or edge[1] == v) return;
            const get = self.map.get(.{ edge[1], edge[0] }).?;
            const w = get.vtx orelse return;
            if (self.inCircumcircle(.{ edge[1], edge[0], w }, v)) {

                    self.deleteTriangle(.{ edge[1], edge[0], w }) catch |err| switch (err) {
                        error.Clobber => return error.Clobber,
                        error.NoTriangle => unreachable,
                    };
                    try self.digCavity(v, .{ edge[0], w });
                    try self.digCavity(v, .{ w, edge[1] });

            } else {
                try self.addTriangle(.{ edge[0], edge[1], v });
            }
        }

        fn findTriangle(self: *const Self, v: usize) [3]usize {
            const keys = self.map.keys();
            var edge: [2]usize, var w = edge: {
                for (0..keys.len) |i| {
                    const edge = keys[keys.len - i - 1];
                    const w = self.map.get(edge).?.vtx orelse continue;
                    break :edge .{ edge, w };
                }
                unreachable;
            };
            while (true) {
                const edges: [3][2]usize = .{
                    .{ edge[1], w },
                    .{ w, edge[0] },
                    .{ edge[0], edge[1] },
                };
                for (edges) |e| {
                    if (!self.inHalfspace(e, v)) {
                        edge = .{ e[1], e[0] };
                        w = self.map.get(edge).?.vtx.?;
                        break;
                    }
                } else return .{ edge[0], edge[1], w };
            }
        }

        fn inCircumcircle(self: *const Self, tri: [3]usize, v: usize) bool {
            if (v == ghost) {
                return false;
            }
            var is_ghost: usize = 3;
            for (tri, 0..) |w, i| {
                if (w == ghost) {
                    is_ghost = i;
                    break;
                }
            }
            if (is_ghost < 3) {
                return switch (is_ghost) {
                    0 => self.inHalfspace(.{ tri[1], tri[2] }, v),
                    1 => return self.inHalfspace(.{ tri[2], tri[0] }, v),
                    2 => return self.inHalfspace(.{ tri[0], tri[1] }, v),
                    else => unreachable,
                };
            }
            const o = math.geometry.orient(2, F, .{
                self.verts.items[tri[0]],
                self.verts.items[tri[1]],
                self.verts.items[tri[2]],
            });
            switch (o) {
                .gt, .eq => return switch (math.geometry.circumsphereContains(2, F, .{
                    self.verts.items[tri[0]],
                    self.verts.items[tri[1]],
                    self.verts.items[tri[2]],
                    self.verts.items[v],
                })) {
                    .gt, .eq => true,
                    .lt => return false,
                },
                .lt => return switch (math.geometry.circumsphereContains(2, F, .{
                    self.verts.items[tri[0]],
                    self.verts.items[tri[1]],
                    self.verts.items[tri[2]],
                    self.verts.items[v],
                })) {
                    .lt => true,
                    .gt, .eq => false,
                },
            }
        }

        fn inHalfspace(self: *const Self, edge: [2]usize, v: usize) bool {
            if (edge[0] == ghost or edge[1] == ghost or v == ghost) return true;
            const o = math.geometry.orient(2, F, .{
                self.verts.items[edge[0]],
                self.verts.items[edge[1]],
                self.verts.items[v],
            });
            return switch (o) {
                .eq => between(self.verts.items[v], self.verts.items[edge[0]], self.verts.items[edge[1]]),
                .lt => false,
                .gt => true,
            };
        }

        fn between(mid: [2]F, a: [2]F, b: [2]F) bool {
            const V = @Vector(2, F);
            const am: V = @as(V, a) - @as(V, mid);
            const bm: V = @as(V, b) - @as(V, mid);
            return switch (math.sign(math.dotProduct(2, F, am, bm))) {
                .eq, .lt => true,
                .gt => false,
            };
        }
    };
}
 
test "planar delaunay" {
    const verts: []const [2]f32 = &.{
        .{ 0, 0 },
        .{ 1, 0 },
        .{ 0, 1 },
        .{ 1, 1 },
    };
    var triangulation = PlanarTriangulation(f32).init(std.testing.allocator);
    defer triangulation.deinit();
    for (verts) |v| {
        try triangulation.addVertex(v);
    }
    const triangles = try triangulation.getNonGhostTriangles(std.testing.allocator);
    defer std.testing.allocator.free(triangles);
    try std.testing.expectEqual(6, triangles.len);
}

test "cool S" {
    const verts: []const [2]f32 = &.{
        .{ -1, -2 },
        .{ 5, 2 },
        .{ 5, 4 },
        .{ 0, 6 },
        .{ -5, 4 },
        .{ -5, 2 },
        .{ -3, 0.25 },
        .{ -1, 2 },
        .{ -1, 4 },
        .{ 1, 4 },
        .{ 1, 2 },
        .{ -5, -2 },
        .{ -5, -4 },
        .{ 0, -6 },
        .{ 5, -4 },
        .{ 5, -2 },
        .{ 3, -0.25 },
        .{ 1, -2 },
        .{ 1, -4 },
        .{ -1, -4 },
    };
    var triangulation = PlanarTriangulation(f32).init(std.testing.allocator);
    defer triangulation.deinit();
    for (verts, 0..) |v, i| {
        if (i < verts.len - 1)
            try triangulation.addEdge(.{ v, verts[i + 1] })
        else
            try triangulation.addEdge(.{ v, verts[0] });
    }
    try triangulation.addEdge(.{verts[0], verts[1]});
    const triangles = try triangulation.getNonGhostTriangles(std.testing.allocator);
    defer std.testing.allocator.free(triangles);
    var oriented_triangles = std.ArrayList(usize).init(std.testing.allocator);
    defer oriented_triangles.deinit();
    var idx: usize = 0;
    while (idx + 3 <= triangles.len) : (idx += 3) {
        const t: [3]usize = triangles[idx..][0..3].*;
        const edges: [3][2]usize = .{
            .{ t[0], t[1] },
            .{ t[1], t[2] },
            .{ t[2], t[0] },
        };
        for (edges) |e| {
            if (e[1] == e[0] + 1 or (e[1] == 0 and e[0] == triangulation.verts.items.len - 1)) {
                try oriented_triangles.appendSlice(&t);
                break;
            }
        }
    }
    // the cool S is a topological disk with 20 vertices, all on the boundary,
    // so because its Euler characteristic is 1, we have 20 - E + F = 1.
    // each triangle has three edges, and 20 edges belong to exactly one triangle
    // so we have 3F = 20 + 2(E - 20) and therefore F = 18
    try std.testing.expectEqual(3 * 18, oriented_triangles.items.len);
}

test "letter 'e'" {
    const points: []const [2]f32 = &.{.{ 2, 5 },
                .{ 0, 5 },
                .{ -2, 5 },
                .{ -3, 4 },
                .{ -4, 3 },
                .{ -4, 0 },
                .{ -4, -3 },
                .{ -3, -4 },
                .{ -2, -5 },
                .{ 0, -5 },
                .{ 2, -5 },
                .{ 3, -4 },
                .{ 3.5, -3 },
                .{ 4, -2.5 },
                .{ 3, -2.5 },
                .{ 2.5, -3 },
                .{ 2, -3.25 },
                .{ 1.5, -3.5 },
                .{ 0, -4 },
                .{ -1.5, -3.5 },
                .{ -2, -3.25 },
                .{ -2.5, -2.5 },
                .{ -3, 0 },
                .{ 4, 0 },
                .{ 4, 0.5 },
                .{ 4, 3 },
                .{ 3, 4 },
                .{ -3, 1 },
                .{ 3, 1 },
                .{ 3, 2 },
                .{ 2.5, 2.5 },
                .{ 1.5, 4 },
                .{ 0, 4 },
                .{ -1.5, 4 },
                .{ -2.5, 2.5 },
                .{ -3, 2 },
            };
            const bezier_triangles: []const struct { [3]usize, bool } = &.{
                     .{ .{ 1, 2, 3 }, true },
                .{ .{ 3, 4, 5 }, true },
                .{ .{ 5, 6, 7 }, true },
                .{ .{ 7, 8, 9 }, true },
                .{ .{ 9, 10, 11 }, true },
                .{ .{ 11, 12, 13 }, true },
                .{ .{ 14, 15, 16 }, false },
                .{ .{ 16, 17, 18 }, false },
                .{ .{ 18, 19, 20 }, false },
                .{ .{ 20, 21, 22 }, false },
                .{ .{ 24, 25, 26 }, true },
                .{ .{ 26, 0, 1 }, true },
                .{ .{ 28, 29, 30 }, false },
                .{ .{ 30, 31, 32 }, false },
                .{ .{ 32, 33, 34 }, false },
                .{ .{ 34, 35, 27 }, false },
            };
           
            var triangulation = PlanarTriangulation(f32).init(std.testing.allocator);
            defer triangulation.deinit();
            for (points) |point| {
                try triangulation.addVertex(point);
            }
            for (bezier_triangles) |tuple| {
                const edges: [3][2]usize = .{
                    .{ tuple[0][0], tuple[0][1] },
                    .{ tuple[0][1], tuple[0][2] },
                    .{ tuple[0][2], tuple[0][0] },
                };
                for (edges) |e| {
                    try triangulation.addEdge(.{ points[e[0]], points[e[1]] });
                }
            }
            const tris = try triangulation.getNonGhostTriangles(std.testing.allocator);
            defer std.testing.allocator.free(tris);
    // var idx: usize = 0;
    // while (idx + 3 <= tris.len) : (idx += 3) {
        // std.debug.print("triangle: {any}\n", .{tris[idx..][0..3]});
    // }
}
