# Replace oneTBB's `libtbbmalloc_proxy` in a PackageCompiler bundle with a
# no-op stub dylib (macOS only).
#
# Why this is necessary
# ---------------------
# FFTW.jl and LinearSolve.jl (a Gradus dependency) depend on MKL_jll, which
# depends on oneTBB_jll, so PackageCompiler bundles the oneTBB libraries into
# `build/share/julia/artifacts/`. On macOS the generated oneTBB_jll wrapper
# eagerly dlopens all three TBB libraries in its `__init__` — including
# `libtbbmalloc_proxy` — and in a PackageCompiler library every bundled
# module's `__init__` runs as soon as XSPEC loads `libgradusxspec.dylib`.
#
# `libtbbmalloc_proxy` exists for exactly one purpose: on load it replaces the
# process-wide macOS default malloc zone with TBB's allocator. Memory that
# XSPEC/Tcl allocated with system malloc *before* our library loads is then
# freed through TBB's allocator *afterwards*, which segfaults, e.g.:
#
#     signal 11: Segmentation fault
#     __TBB_malloc_safer_msize at .../libtbbmalloc.2.17.dylib
#     (crash report: find_zone_and_free -> __TBB_malloc_safer_msize inside
#      xs_execute_script while xspec sources its RC file)
#
# Nothing in the process ever calls TBB (MKL does not even exist on Apple
# Silicon), so the proxy is pure collateral damage. There is no runtime
# kill-switch on macOS (TBB_MALLOC_DISABLE_REPLACEMENT is Windows-only) and
# MKL_jll cannot be dropped from the manifest, so we swap the bundled proxy
# for an empty dylib with the same install name. oneTBB_jll's `__init__`
# still dlopens it successfully, but no malloc-zone replacement happens.
#
# Linux is unaffected: the JLL dlopens without RTLD_GLOBAL, so the proxy's
# malloc symbols never interpose there.
#
# Run automatically from `src/build_lib.jl` after `create_library`, or
# standalone against an existing bundle:
#
#     julia scripts/stub_tbbmalloc_proxy.jl [build_dir]

function stub_tbbmalloc_proxy!(build_dir::AbstractString = "build")
    if !Sys.isapple()
        println("stub_tbbmalloc_proxy: not macOS; nothing to do")
        return 0
    end
    artifacts = joinpath(build_dir, "share", "julia", "artifacts")
    if !isdir(artifacts)
        println("stub_tbbmalloc_proxy: no artifacts directory at $artifacts")
        return 0
    end

    stub_src = tempname() * ".c"
    write(stub_src, "/* intentionally empty: stub for libtbbmalloc_proxy */\n")

    n = 0
    for (root, _dirs, files) in walkdir(artifacts)
        for f in files
            startswith(f, "libtbbmalloc_proxy") && endswith(f, ".dylib") || continue
            path = joinpath(root, f)
            islink(path) && continue  # symlinks resolve to the stubbed file
            stub = path * ".stub.dylib"
            run(`cc -dynamiclib -o $stub $stub_src -install_name @rpath/libtbbmalloc_proxy.2.dylib`)
            # Bundled artifact files are read-only; replace the directory entry.
            rm(path)
            mv(stub, path)
            chmod(path, 0o555)
            println("stub_tbbmalloc_proxy: replaced $path")
            n += 1
        end
    end
    rm(stub_src; force = true)
    if n == 0
        println("stub_tbbmalloc_proxy: no libtbbmalloc_proxy dylibs found (ok if oneTBB is not bundled)")
    end
    return n
end

if abspath(PROGRAM_FILE) == @__FILE__
    stub_tbbmalloc_proxy!(isempty(ARGS) ? "build" : ARGS[1])
end
