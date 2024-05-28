const math = @import("math.zig");
const std = @import("std");
pub const delaunay = @import("delaunay.zig");

/// determines whether two points are approximately equal
pub fn eqPts(comptime n: usize, comptime F: type, a: [n]F, b: [n]F) bool {
    var ret = true;
    // is unrolling a good idea?
    inline for (0..n) |i| {
        ret = ret and math.eq(a[i], b[i]);
    }
    return ret;
}

fn affineLift(comptime m: usize, comptime n: usize, comptime F: type, pts: [m][n]F) [m + 1][n]F {
    var ret: [m + 1][n]F = .{.{1} ** n} ** (m + 1);
    @memcpy(ret[0..m], &pts);
    return ret;
}

/// determines the orientation of the affine N-simplex spanned by the N+1 points.
pub fn orient(comptime n: usize, comptime F: type, pts: [n + 1][n]F) std.math.Order {
    const mat = affineLift(n, n + 1, F, math.transpose(n + 1, n, F, pts));
    return math.sign(math.determinant(n + 1, F, mat));
}

/// determines whether the circumsphere of the first  N+1 points in N-space
/// contains the last point
/// the first N+1 points are assumed positively oriented.
pub fn circumsphereContains(comptime n: usize, comptime F: type, points: [n + 2][n]F) std.math.Order {
    var mat: [n + 2][n + 2]F = .{.{1} ** (n + 2)} ** (n + 2);
    @memcpy(mat[0..n], &math.transpose(n + 2, n, F, points));
    inline for (0..n + 2) |i| {
        mat[n][i] = math.dotProduct(n, F, points[i], points[i]);
    }
    const det = math.determinant(n + 2, F, mat);
    return math.sign(det);
}

/// determines whether the given affine subspaces intersect
/// returns a list of points affinely generating the intersection if so
pub fn affineIntersection(
    comptime n: usize,
    comptime a: usize,
    comptime b: usize,
    comptime F: type,
    A: [a][n]F,
    B: [b][n]F,
) struct { [a + b][n]F, usize } {
    const ABt = math.transpose(a + b, n, F, A ++ B);
    var matrix: [n + 1][a + b]F = affineLift(n, a + b, F, ABt);
    inline for (0..b) |i| {
        inline for (0..n + 1) |j| {
            matrix[j][a + i] *= -1;
        }
    }
    var ret: [a + b][n]F = .{.{0} ** n} ** (a + b);
    const null_space, const i = math.nullSpace(n + 1, a + b, F, matrix);
    var curr: usize = 0;
    for (0..i) |j| {
        const w = w: {
            var w: F = 0;
            inline for (0..a + b) |k| {
                w += null_space[j][k];
            }
            break :w w;
        };
        switch (math.sign(w)) {
            .eq => continue,
            .gt, .lt => {
                inline for (0..n) |k| {
                    ret[curr][k] = math.dotProduct(a + b, F, ABt[k], null_space[j]) / w;
                }
                curr += 1;
            },
        }
    }
    return .{ ret, curr };
}

pub fn triangulatePlanarQuad(comptime F: type, verts: [4][2]F) [6]usize {
    const edges: [6]struct { [2]usize, [2]usize } = .{
        .{ .{ 0, 1 }, .{ 2, 3 } },
        .{ .{ 0, 2 }, .{ 1, 3 } },
        .{ .{ 0, 3 }, .{ 1, 2 } },
        .{ .{ 1, 2 }, .{ 0, 3 } },
        .{ .{ 1, 3 }, .{ 0, 2 } },
        .{ .{ 2, 3 }, .{ 0, 1 } },
    };
    for (edges) |diag| {
        const a = orient(2, F, .{
            verts[diag[0][0]],
            verts[diag[0][1]],
            verts[diag[1][0]],
        });
        const b = orient(2, F, .{
            verts[diag[0][0]],
            verts[diag[0][1]],
            verts[diag[1][1]],
        });
            if (a != b) {
                return .{ diag[1][0], diag[0][0], diag[0][1], diag[0][0], diag[0][1], diag[1][1], };
        }
    }
    return .{
        0, 1, 2, 2, 1, 3 
    };
}

test "intersection" {
    const A: [2][2]f32 = .{
        .{ 1, 3 },
        .{ 3, 1 },
    };
    const B: [2][2]f32 = .{
        .{ 3, 3 },
        .{ 1, 1 },
    };
    const C: [2][2]f32 = .{
        .{ 5, 3 },
        .{ 3, 1 },
    };
    const pts, const i = affineIntersection(2, 2, 2, f32, A, B);
    try std.testing.expectEqual(1, i);
    for (pts[0]) |pt| {
        try std.testing.expectApproxEqAbs(2, pt, math.eps);
    }
    _, const j = affineIntersection(2, 2, 2, f32, B, C);
    try std.testing.expectEqual(0, j);
}

test "orientation" {
    const verts: [3][2]f32 = .{
        .{ 0, 0 },
        .{ 1, 0 },
        .{ 0, 1 },
    };
    try std.testing.expectEqual(.gt, orient(2, f32, verts));

    const verts2: [4][3]f32 = .{
        .{ 1, 0, 0 },
        .{ 0, 1, 0 },
        .{ 0, 0, 0 },
        .{ 0, 0, 1 },
    };
    try std.testing.expectEqual(.lt, orient(3, f32, verts2));
}

test "circumsphere" {
    const verts: [4][2]f32 = .{
        .{ 0, 0 },
        .{ 1, 0 },
        .{ 0, 2 },
        .{ 1, 2 },
    };
    try std.testing.expectEqual(.eq, circumsphereContains(2, f32, verts));

    const verts2: [4][2]f32 = .{
        .{ 2, -2 },
        .{ 0, 4 },
        .{ -2, -2 },
        .{ 0, 0 },
    };
    try std.testing.expectEqual(.gt, circumsphereContains(2, f32, verts2));
}
