const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("math", .{
        .target = target,
        .optimize = optimize,
        .root_source_file = .{ .path = "src/math.zig" },
    });

    const tests = b.addTest(.{
        .root_source_file = .{ .path = "src/math.zig" },
        .target = target,
        .optimize = optimize,
    });
    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "run tests");
    test_step.dependOn(&run_tests.step);
}
