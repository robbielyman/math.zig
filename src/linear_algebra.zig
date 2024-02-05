const std = @import("std");
const math = @import("math.zig");

/// finds the norm (length) of a vector
pub fn length(comptime n: usize, comptime F: type, a: [n]F) F {
    return @sqrt(dotProduct(n, F, a, a));
}

/// performs the dot product of two vectors
pub fn dotProduct(comptime n: usize, comptime F: type, a: [n]F, b: [n]F) F {
    const V = @Vector(n, F);
    return @reduce(.Add, @as(V, a) * @as(V, b));
}

/// transposes a rectangular matrix.
pub fn transpose(comptime m: usize, comptime n: usize, comptime T: type, matrix: [m][n]T) [n][m]T {
    var ret: [n][m]T = undefined;
    // is unrolling a good idea?
    inline for (0..n) |i| {
        inline for (0..m) |j| {
            ret[i][j] = matrix[j][i];
        }
    }
    return ret;
}

/// `matrix` is interpreted as an array of _row_ vectors.
/// F should be a floating point type: division is performed.
pub fn determinant(comptime n: usize, comptime F: type, matrix: [n][n]F) F {
    const V = @Vector(n, F);
    var mat: [n]V = undefined;
    for (0..n) |i| {
        mat[i] = matrix[i];
    }
    var det: F = 1.0;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        // the next "pivot" could be any nonzero entry,
        // but numerically we'd prefer to divide by the largest number.
        const j = find_pivot: {
            var pivot: usize = n + 1;
            var max: F = 0;
            for (i..n) |j| {
                if (math.order(@abs(mat[j][i]), max) == .gt) {
                    max = @abs(mat[j][i]);
                    pivot = j;
                }
            }
            // if a whole column is zeo, the determinant is zero.
            if (math.sign(max) == .gt) break :find_pivot pivot else return 0;
        };
        // swap rows
        if (j != i) {
            const v: V = mat[j];
            mat[j] = mat[i];
            mat[i] = v;
            // swapping multiplies the determinant by -1
            det = -det;
        }
        // rescale the row so the pivot entry is 1.
        // scaling a row multiplies the determinant by the reciprocal of the scaling factor
        det *= mat[i][i];
        mat[i] /= @splat(mat[i][i]);
        // zero out the rest of the column by performing row operations.
        // numerically, I suppose we could ignore the columns that are already small?
        for (i + 1..n) |k| {
            const v: V = @splat(mat[k][i]);
            mat[k] -= v * mat[i];
        }
    }
    return det;
}

/// `mat` is interpreted as an array of _row_ vectors.
/// F should be a floating point type: division is performed.
pub fn rowReduce(comptime n: usize, comptime F: type, matrix: [n][n]F) [n][n]F {
    const V = @Vector(n, F);
    var mat: [n]V = undefined;
    for (0..n) |i| {
        mat[i] = matrix[i];
    }
    var i: usize = 0;
    while (i < n) : (i += 1) {
        // the next "pivot" could be any nonzero entry,
        // but numerically we'd prefer to divide by the largest number;
        const j: usize = find_pivot: {
            var pivot: usize = 0;
            var max: F = -1.0;
            for (i..n) |j| {
                if (math.order(@abs(mat[j][i]), max) == .gt) {
                    max = @abs(mat[j][i]);
                    pivot = j;
                }
            }
            if (math.sign(max) == .gt) break :find_pivot pivot else continue;
        };
        // swap rows
        if (j != i) {
            const v: V = mat[j];
            mat[j] = mat[i];
            mat[i] = v;
        }
        // rescale the row so the pivot entry is 1.
        mat[i] /= @splat(mat[i][i]);
        // zero out the rest of the column by performing row operations.
        // numerically, I suppose we could ignore the columns that are already small?
        for (i + 1..n) |k| {
            const v: V = @splat(mat[k][i]);
            mat[k] -= v * mat[i];
        }
    }
    var ret: [n][n]F = undefined;
    for (0..n) |j| {
        ret[j] = mat[j];
    }
    return ret;
}

test "determinants" {
    const testing = std.testing;
    const A: [2][2]f32 = .{
        .{ 1, 1 },
        .{ 2, 2 },
    };
    try testing.expectApproxEqAbs(0, determinant(2, f32, A), math.eps);
    const B: [3][3]f32 = .{
        .{ 1, 2, 3 },
        .{ 0, 1, 2 },
        .{ 0, 0, 1 },
    };
    try testing.expectApproxEqAbs(1, determinant(3, f32, B), math.eps);
    const C: [2][2]f32 = .{
        .{ math.eps / 2.0, 0.5 },
        .{ 2, 1 },
    };
    try testing.expectApproxEqAbs(-1, determinant(2, f32, C), math.eps);
}

test "row reduce" {
    const testing = std.testing;
    const A: [2][2]f32 = .{
        .{ 1, 1 },
        .{ 2, 2 },
    };
    const red_A = rowReduce(2, f32, A);
    for (0..2) |i| {
        try testing.expectApproxEqAbs(0, red_A[1][i], math.eps);
    }
    const B: [3][3]f32 = .{
        .{ 1, 2, 3 },
        .{ 0, 1, 2 },
        .{ 0, 0, 1 },
    };
    const red_B = rowReduce(3, f32, B);
    for (0..3) |i| {
        try std.testing.expectEqualSlices(f32, &B[i], &red_B[i]);
    }
    const C: [2][2]f32 = .{
        .{ math.eps / 2.0, 0.5 },
        .{ 2, 1 },
    };
    const red_C = rowReduce(2, f32, C);
    const expected: [2][2]f32 = .{
        .{ 1, 0.5 },
        .{ 0, 1 },
    };
    for (red_C, expected) |row, ex| {
        for (row, ex) |actual, e| {
            try testing.expectApproxEqAbs(e, actual, math.eps);
        }
    }
}
