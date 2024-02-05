const std = @import("std");
const math = @import("math.zig");

// the triangulation is a list of edges,
// with each edge pointing to the two triangles containing it
// we work with a "point at infinity", represented by index n,
// so every "real" edge belongs to exactly two triangles,
const DelaunayMap = std.AutoArrayHashMap([2]usize, TriangleData);

const TriangleData = struct {
    v: ?usize,
    flippable: bool = true,
};

pub const Err = error{ Degenerate, Clobber, NoTriangle } || std.mem.Allocator.Error;

fn PlanarTriangulation(comptime F: type) type {
    return struct {
        const Self = @This();
        map: DelaunayMap,
        verts: []const [2]F,

        fn addInitialTriangle(self: *Self, tri: [3]usize) Err!void {
            const n = self.verts.len;
            switch (math.geometry.sphericalOrder(2, F, .{
                self.verts[tri[0]],
                self.verts[tri[1]],
                self.verts[tri[2]],
            })) {
                .gt => {
                    try self.addTriangle(.{ tri[0], tri[1] }, tri[2]);
                    try self.addTriangle(.{ tri[1], tri[0] }, n);
                    try self.addTriangle(.{ tri[0], tri[2] }, n);
                    try self.addTriangle(.{ tri[2], tri[1] }, n);
                },
                .lt => {
                    try self.addTriangle(.{ tri[1], tri[0] }, tri[2]);
                    try self.addTriangle(.{ tri[0], tri[1] }, n);
                    try self.addTriangle(.{ tri[1], tri[2] }, n);
                    try self.addTriangle(.{ tri[2], tri[0] }, n);
                },
                .eq => {
                    const eq = math.geometry.eqPts(2, F, self.verts[tri[0]], self.verts[tri[1]]) or
                        math.geometry.eqPts(2, F, self.verts[tri[0]], self.verts[tri[2]]) or
                        math.geometry.eqPts(2, F, self.verts[tri[1]], self.verts[tri[2]]);
                    if (eq) return error.Degenerate;
                    if (try self.between(tri[1], tri[0], tri[2])) {
                        try self.addTriangle(.{ tri[0], tri[1] }, tri[2]);
                        try self.addTriangle(.{ tri[1], tri[0] }, n);
                        try self.addTriangle(.{ tri[0], tri[2] }, n);
                        try self.addTriangle(.{ tri[2], tri[1] }, n);
                    } else {
                        try self.addTriangle(.{ tri[1], tri[0] }, tri[2]);
                        try self.addTriangle(.{ tri[0], tri[1] }, n);
                        try self.addTriangle(.{ tri[1], tri[2] }, n);
                        try self.addTriangle(.{ tri[2], tri[0] }, n);
                    }
                },
            }
        }

        fn between(self: *const Self, mid: usize, a: usize, b: usize) Err!bool {
            const V = @Vector(2, F);
            const v = @as(V, self.verts[a]) - @as(V, self.verts[mid]);
            const w = @as(V, self.verts[b]) - @as(V, self.verts[mid]);
            return switch (math.sign(math.dotProduct(2, F, v, w))) {
                .eq => error.Degenerate,
                .gt => true,
                .lt => false,
            };
        }

        /// edge is assumed to be positively oriented with respect to v
        fn addTriangle(self: *Self, edge: [2]usize, v: usize) Err!void {
            const edges: [3][2]usize = .{
                .{ edge[0], edge[1] },
                .{ edge[1], v },
                .{ v, edge[0] },
            };
            const other: [3]usize = .{ v, edge[0], edge[1] };
            for (edges) |e| {
                const res = try self.map.getOrPut(e);
                if (res.found_existing) {
                    if (res.value_ptr.v != null) return error.Clobber;
                } else {
                    res.value_ptr.* = .{ .v = null };
                }
            }
            for (edges, other) |e, w| {
                const ptr = self.map.getPtr(e).?;
                ptr.v = w;
            }
        }

        fn deleteTriangle(self: *Self, tri: [3]usize) Err!void {
            const edges: [3][2]usize = .{ .{ tri[1], tri[2] }, .{ tri[2], tri[0] }, .{ tri[0], tri[1] } };
            for (edges, 0..) |edge, i| {
                const val = self.map.getPtr(edge) orelse return error.NoTriangle;
                if (val.v) |w| if (w != tri[i]) return error.Clobber;
                val.v = null;
            }
        }

        fn inCircle(self: *const Self, tri: [3]usize, v: usize) Err!bool {
            const n = self.verts.len;
            var ghost: usize = 3;
            for (tri, 0..) |w, i| {
                if (w == n) ghost = i;
            }
            if (ghost < 3) {
                const edge: [2]usize = switch (ghost) {
                    0 => .{ tri[1], tri[2] },
                    1 => .{ tri[2], tri[0] },
                    2 => .{ tri[0], tri[1] },
                    else => unreachable,
                };
                return switch (math.geometry.sphericalOrder(2, F, .{ self.verts[edge[0]], self.verts[edge[1]], self.verts[v] })) {
                    .gt => true,
                    .lt => false,
                    .eq => try self.between(v, edge[0], edge[1]),
                };
            }
            return switch (math.geometry.circumsphereContains(2, F, .{
                self.verts[tri[0]],
                self.verts[tri[1]],
                self.verts[tri[2]],
                self.verts[v],
            })) {
                .gt, .eq => true,
                .lt => false,
            };
        }

        fn findTriangle(self: *const Self, first_edge: [2]usize, v: usize) Err![3]usize {
            const n = self.verts.len;
            var edge = first_edge;
            while (true) {
                if (edge[0] == n or edge[1] == n) {
                    const w = self.map.get(edge).?.v.?;
                    edge = .{ edge[1], w };
                    continue;
                }
                switch (math.geometry.sphericalOrder(2, F, .{
                    self.verts[edge[0]], self.verts[edge[1]], self.verts[v],
                })) {
                    .gt => {
                        const w = self.map.get(edge).?.v.?;
                        const tri: [3]usize = .{ edge[0], edge[1], w };
                        if (try self.inCircle(tri, v)) return tri;
                        edge = .{ edge[1], w };
                    },
                    .lt => {
                        const w1 = self.map.get(edge).?.v.?;
                        const data = self.map.get(.{ edge[1], edge[0] });
                        if (data) |d| {
                            const w = d.v.?;
                            const tmp = .{ edge[1], edge[0] };
                            edge = tmp;
                            const tri: [3]usize = .{ edge[1], edge[0], w };
                            if (try self.inCircle(tri, v)) return tri;
                            edge = .{ edge[1], w };
                        } else edge = .{ edge[1], w1 };
                    },
                    .eq => {
                        const w = self.map.get(edge).?.v.?;
                        if (try self.between(v, edge[0], edge[1])) {
                            return .{ edge[0], edge[1], w };
                        }
                        edge = .{ edge[1], w };
                    },
                }
            }
        }

        fn digCavity(self: *Self, edge: [2]usize, v: usize) Err!void {
            const data = self.map.getPtr(.{ edge[1], edge[0] }).?;
            const w = data.v orelse return;
            if (try self.inCircle(.{ edge[1], edge[0], w }, v)) {
                try self.deleteTriangle(.{ edge[1], edge[0], w });
                try self.digCavity(.{ edge[0], w }, v);
                try self.digCavity(.{ w, edge[1] }, v);
            } else {
                try self.addTriangle(edge, v);
            }
        }

        fn deinit(self: *Self) void {
            self.map.deinit();
        }
    };
}

/// returns a list of indices into verts;
/// the length of this list will be a multiple of three (unless verts is too small)
/// and each group of three determines a triangle of the triangulation
pub fn planarDelaunay(comptime F: type, allocator: std.mem.Allocator, verts: []const [2]F) Err![]usize {
    const n = verts.len;
    const indices = try biasedRandomOrder(n, allocator);
    defer allocator.free(indices);
    var triangles = std.ArrayList(usize).init(allocator);
    defer triangles.deinit();
    // there is at most one triangle, so just return that triangle
    if (n < 4) {
        try triangles.appendSlice(indices);
        return try triangles.toOwnedSlice();
    }

    var triangulation: PlanarTriangulation(F) = .{
        .map = DelaunayMap.init(allocator),
        .verts = verts,
    };
    defer triangulation.deinit();
    try triangulation.addInitialTriangle(.{ indices[0], indices[1], indices[2] });

    for (indices[3..], 3..) |v, i| {
        const edge = edge: {
            for (triangulation.map.keys()) |k| {
                if (k[0] == indices[i - 1] or k[1] == indices[i - 1]) break :edge k;
            }
            unreachable;
        };
        const tri = try triangulation.findTriangle(edge, v);
        try triangulation.deleteTriangle(tri);
        try triangulation.digCavity(.{ tri[0], tri[1] }, v);
        try triangulation.digCavity(.{ tri[1], tri[2] }, v);
        try triangulation.digCavity(.{ tri[2], tri[0] }, v);
    }

    for (triangulation.map.keys()) |edge| {
        if (triangulation.map.getPtr(edge).?.v) |v| {
            if (v != n and edge[0] != n and edge[1] != n) {
                try triangles.appendSlice(&.{ edge[0], edge[1], v });
            }
            try triangulation.deleteTriangle(.{ edge[0], edge[1], v });
        }
    }
    return try triangles.toOwnedSlice();
}

fn biasedRandomOrder(n: usize, allocator: std.mem.Allocator) std.mem.Allocator.Error![]usize {
    var rand = std.rand.DefaultPrng.init(@intCast(@max(0, std.time.nanoTimestamp())));
    const r = rand.random();
    const list = try allocator.alloc(usize, n);
    var idx: usize = 0;
    while (idx < n) {
        for (0..n) |i| {
            for (0..idx) |j| {
                if (list[j] == i) break;
            } else {
                const flip = r.int(u1);
                switch (flip) {
                    0 => {},
                    1 => {
                        list[idx] = i;
                        idx += 1;
                    },
                }
            }
        }
    }
    return list;
}

test "planar delaunay" {
    const verts: []const [2]f32 = &.{
        .{ 0, 0 },
        .{ 1, 0 },
        .{ 0, 1 },
        .{ 1, 1 },
    };
    const triangulation = try planarDelaunay(f32, std.testing.allocator, verts);
    defer std.testing.allocator.free(triangulation);
    try std.testing.expectEqual(6, triangulation.len);
}
