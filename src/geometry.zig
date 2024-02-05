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

/// determines the "spherical order" of the N+1 points in N-space
/// that is, thinking of N-space projectively
/// do they, with the origin, form a positively oriented simplex?
pub fn sphericalOrder(comptime n: usize, comptime F: type, pts: [n + 1][n]F) std.math.Order {
    var mat: [n + 1][n + 1]F = .{.{1} ** (n + 1)} ** (n + 1);
    @memcpy(mat[1 .. n + 1], &math.transpose(n + 1, n, F, pts));
    return math.sign(math.determinant(n + 1, F, mat));
}

/// determines whether the circumsphere of the first  N+1 points in N-space
/// contains the last point
pub fn circumsphereContains(comptime n: usize, comptime F: type, points: [n + 2][n]F) std.math.Order {
    const orientation = sphericalOrder(n, F, points[0 .. n + 1].*);
    var mat: [n + 2][n + 2]F = .{.{1} ** (n + 2)} ** (n + 2);
    @memcpy(mat[0..n], &math.transpose(n + 2, n, F, points));
    inline for (0..n + 2) |i| {
        mat[n][i] = math.dotProduct(n, F, points[i], points[i]);
    }
    const inside = math.sign(math.determinant(n + 2, F, mat));
    return switch (orientation) {
        .gt => inside,
        .lt => math.reverse(inside),
        .eq => .eq,
    };
}

test "spherical order" {
    const verts: [3][2]f32 = .{
        .{ 0, 0 },
        .{ 1, 0 },
        .{ 0, 1 },
    };
    try std.testing.expectEqual(.gt, sphericalOrder(2, f32, verts));

    const verts2: [4][3]f32 = .{
        .{ 1, 0, 0 },
        .{ 0, 0, 0 },
        .{ 0, 1, 0 },
        .{ 0, 0, 1 },
    };
    try std.testing.expectEqual(.lt, sphericalOrder(3, f32, verts2));
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
        .{ -2, -2 },
        .{ 0, 4 },
        .{ 0, 0 },
    };
    try std.testing.expectEqual(.gt, circumsphereContains(2, f32, verts2));
}
