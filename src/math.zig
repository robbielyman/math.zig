const std = @import("std");
const linear_algebra = @import("linear_algebra.zig");
pub const geometry = @import("geometry.zig");
pub const cubic = @import("cubic.zig");

pub const eps: comptime_float = 0.000001;

pub const determinant = linear_algebra.determinant;
pub const rowReduce = linear_algebra.rowReduce;
pub const transpose = linear_algebra.transpose;
pub const dotProduct = linear_algebra.dotProduct;
pub const crossProduct = linear_algebra.crossProduct;
pub const nullSpace = linear_algebra.nullSpace;
pub const length = linear_algebra.length;
pub const multiply = linear_algebra.multiply;
pub const dimensionOfSpan = linear_algebra.dimensionOfSpan;
pub const identityMatrix = linear_algebra.identity;

/// fuzzy ordering; differences of size `eps` are treated as equality.
/// `a` and `b` should be floating point types
pub fn order(a: anytype, b: anytype) std.math.Order {
    if (@abs(a - b) <= eps) return .eq;
    return if (a - b > eps) .gt else .lt;
}

/// reverses order
pub fn reverse(o: std.math.Order) std.math.Order {
    return switch (o) {
        .gt => .lt,
        .eq => .eq,
        .lt => .gt,
    };
}

/// fuzzy sign
pub fn sign(a: anytype) std.math.Order {
    if (@abs(a) <= eps) return .eq;
    return if (a > eps) .gt else .lt;
}

/// fuzzy equality check; differences of size `eps` are treated as equality.
/// `a` and `b` should be floating point types
pub fn eq(a: anytype, b: anytype) bool {
    return @abs(a - b) <= eps;
}

test "RefAllDecls" {
    _ = @This();
    try std.testing.expect(eps >= 0);
    std.testing.refAllDeclsRecursive(@This());
}
