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

/// performs the three-dimensional cross product
pub fn crossProduct(comptime F:type, a: [3]F, b: [3]F) [3]F {
    return .{ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],  a[0] * b[1] - b[0] * a[1] };
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

/// determines the dimension of the space the given collection of N-dimensional vectors spans
pub fn dimensionOfSpan(comptime m: usize, comptime n: usize, comptime F: type, matrix: [m][n]F) usize {
    var mat = matrix;
    var ret: usize = 0;
    for (0..m) |i| {
        // the next "pivot" could be any nonzero entry,
        // but numerically we'd prefer to divide by the largest number;
        const j: usize = find_pivot: {
            var pivot: usize = 0;
            var max: F = -1.0;
            for (i..n) |j| {
                if (math.order(@abs(mat[i][j]), max) == .gt) {
                    max = @abs(mat[i][j]);
                    pivot = j;
                }
            }
            if (math.sign(max) != .gt) continue else break :find_pivot pivot;
        };
        ret += 1;
        // swap rows
        if (j != i) {
            inline for (0..m) |k| {
                const v = mat[k][j];
                mat[k][j] = mat[k][i];
                mat[k][i] = v;
            }
        }
        // rescale the row
        const v = mat[i][i];
        inline for (0..m) |k| {
            mat[k][i] /= v;
        }
        // zero out the rest of the column below us
        // we only need to find echelon form
        for (i + 1..n) |k| {
            const w = mat[i][k];
            for (0..m) |l| {
                mat[l][k] -= w * mat[l][i];
            }
        }
    }
    return ret;
}

/// `mat` is interpreted as an array of _row_ vectors.
/// F should be a floating point type: division is performed.
pub fn rowReduce(comptime m: usize, comptime n: usize, comptime F: type, matrix: [m][n]F) [m][n]F {
    const V = @Vector(n, F);
    var mat: [m]V = undefined;
    for (0..m) |i| {
        mat[i] = matrix[i];
    }
    var i: usize = 0;
    while (i < m) : (i += 1) {
        // the next "pivot" could be any nonzero entry,
        // but numerically we'd prefer to divide by the largest number;
        const j: usize = find_pivot: {
            var pivot: usize = 0;
            var max: F = -1.0;
            for (i..m) |j| {
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
        for (0..m) |k| {
            if (k == i) continue;
            const v: V = @splat(mat[k][i]);
            mat[k] -= v * mat[i];
        }
    }
    var ret: [m][n]F = undefined;
    for (0..m) |k| {
        ret[k] = mat[k];
    }
    return ret;
}

/// computes a basis for the null space of the given matrix
/// returns an NxN matrix R and a number k;
/// a basis for the null space is given by R[0..k].
pub fn nullSpace(
    comptime m: usize,
    comptime n: usize,
    comptime F: type,
    matrix: [m][n]F,
) struct { [n][n]F, usize } {
    var mat: [m + n][n]F = matrix ++ .{.{0} ** n} ** n;
    for (0..n) |i| {
        mat[m + i][i] = 1;
    }
    var i: usize = 0;
    while (i < m) : (i += 1) {
        // this is _column_ reduction
        // the next "pivot" could be any nonzero entry,
        // but numerically we'd prefer to divide by the largest number;
        const j: usize = find_pivot: {
            var pivot: usize = 0;
            var max: F = -1.0;
            for (i..n) |j| {
                if (math.order(@abs(mat[i][j]), max) == .gt) {
                    max = @abs(mat[i][j]);
                    pivot = j;
                }
            }
            if (math.sign(max) == .gt) break :find_pivot pivot else continue;
        };
        // swap columns
        if (j != i) {
            inline for (0..m + n) |k| {
                const v = mat[k][j];
                mat[k][j] = mat[k][i];
                mat[k][i] = v;
            }
        }
        // rescale so the pivot entry is 1
        const div = mat[i][i];
        inline for (0..m + n) |k| {
            mat[k][i] /= div;
        }
        // zero out the row by performing column operations
        for (0..n) |l| {
            if (l == i) continue;
            const v = mat[i][l];
            inline for (0..m + n) |k| {
                mat[k][l] -= v * mat[k][i];
            }
        }
    }
    var ret: [n][n]F = .{.{0} ** n} ** n;
    i = 0;
    for (0..n) |j| {
        for (0..m) |k| {
            if (math.sign(mat[k][j]) != .eq) break;
        } else {
            inline for (0..n) |k| {
                ret[i][k] = mat[m + k][j];
            }
            i += 1;
        }
    }
    return .{ ret, i };
}

/// performs matrix multiplication extremely naively.
pub fn multiply(
    comptime m: usize,
    comptime k: usize,
    comptime n: usize,
    comptime F: type,
    A: [m][k]F,
    B: [k][n]F,
) [m][n]F {
    var ret: [m][n]F = undefined;
    for (0..m) |i| {
        for (0..n) |j| {
            const v: [k]F = v: {
                var v: [k]F = undefined;
                inline for (0..k) |l| {
                    v[l] = B[l][j];
                }
                break :v v;
            };
            ret[i][j] = dotProduct(k, F, A[i], v);
        }
    }
    return ret;
}

pub fn identity(comptime n: usize, comptime F: type) [n][n]F {
    var ret: [n][n]F = .{.{0} ** n} ** n;
    inline for (0..n) |i| {
        ret[i][i] = 1;
    }
    return ret;
}

test "linear independence" {
    const testing = std.testing;
    const A: [4][3]f32 = .{
        .{ 3, 1, 1 },
        .{ 1, 3, 1 },
        .{ 3, 3, 1 },
        .{ 1, 1, 1 },
    };
    try testing.expectEqual(3, dimensionOfSpan(4, 3, f32, A));
}

test "null space" {
    const testing = std.testing;
    const A: [3][4]f32 = .{
        .{ 3, 1, -3, -1 },
        .{ 1, 3, -3, -1 },
        .{ 1, 1, -1, -1 },
    };
    const null_space, const len = nullSpace(3, 4, f32, A);
    try testing.expectEqual(1, len);
    for (null_space[0]) |val| {
        try testing.expectApproxEqAbs(1, val, math.eps);
    }
    const B: [2][3]f32 = .{
        .{ 1, 3, 9 },
        .{ 0, 0, 0 },
    };
    const b_null_space, const b_len = nullSpace(2, 3, f32, B);
    try testing.expectEqual(2, b_len);
    const expected: [2][3]f32 = .{ .{ 0, 1, -1.0/3.0 }, .{ 1, 0, -1.0/9.0 } };
    for (b_null_space[0..b_len], expected) |got, expect| {
        for (got, expect) |a, e| {
            try testing.expectApproxEqAbs(e, a, math.eps);
        }
    }
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
    const red_A = rowReduce(2, 2, f32, A);
    for (0..2) |i| {
        try testing.expectApproxEqAbs(0, red_A[1][i], math.eps);
    }
    const B: [3][3]f32 = .{
        .{ 1, 2, 3 },
        .{ 0, 1, 2 },
        .{ 0, 0, 1 },
    };
    const id: [3][3]f32 = .{
        .{ 1, 0, 0 },
        .{ 0, 1, 0 },
        .{ 0, 0, 1 },
    };
    const red_B = rowReduce(3, 3, f32, B);
    for (0..3) |i| {
        try std.testing.expectEqualSlices(f32, &id[i], &red_B[i]);
    }
    const C: [2][2]f32 = .{
        .{ math.eps / 2.0, 0.5 },
        .{ 2, 1 },
    };
    const red_C = rowReduce(2, 2, f32, C);
    const expected: [2][2]f32 = .{
        .{ 1, 0 },
        .{ 0, 1 },
    };
    for (red_C, expected) |row, ex| {
        for (row, ex) |actual, e| {
            try testing.expectApproxEqAbs(e, actual, math.eps);
        }
    }
    const D: [3][4]f32 = .{
        .{ 3, 1, -3, -1 },
        .{ 1, 3, -3, -1 },
        .{ 1, 1, -1, -1 },
    };
    const exp_D: [3][4]f32 = .{
        .{ 1, 0, 0, -1 },
        .{ 0, 1, 0, -1 },
        .{ 0, 0, 1, -1 },
    };
    const red_D = rowReduce(3, 4, f32, D);
    for (red_D, exp_D) |row, ex| {
        for (row, ex) |actual, e| {
            try testing.expectApproxEqAbs(e, actual, math.eps);
        }
    }
}
