# Copyright (c) 2024 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MIT

function decompress_binary(infile, outfile; chunksize=4096)
    open(infile) do input
        open(outfile, "w") do output
            stream = XzDecompressorStream(input)
            while !eof(stream)
                write(output, read(stream, chunksize))
            end
        end
    end
end

function compress_binary(infile, outfile; chunksize=4096)
    open(infile, "r") do input
        open(outfile, "w") do output
            stream = XzCompressorStream(output)
            while !eof(input)
                data = read(input, chunksize)
                write(stream, data)
            end
            close(stream)  # important to flush and close compressor stream
        end
    end
end

function filecmp(path1::AbstractString, path2::AbstractString)
    stat1, stat2 = stat(path1), stat(path2)
    if !(isfile(stat1) && isfile(stat2)) || filesize(stat1) != filesize(stat2)
        return false
    end
    open(path1, "r") do file1
        open(path2, "r") do file2
            buf1 = Vector{UInt8}(undef, 32768)
            buf2 = similar(buf1)
            while !eof(file1) && !eof(file2)
                n1 = readbytes!(file1, buf1)
                n2 = readbytes!(file2, buf2)
                if n1 != n2 || buf1[1:n1] != buf2[1:n2]
                    return false
                end
            end
            return eof(file1) && eof(file2)
        end
    end
end

function create_default_models(; prn=true)
    function create_model(name; segments=3)
        set = Settings("system.yaml")
        set.segments = segments
        set.physical_model = name
        s = SymbolicAWEModel(set)
        time = @elapsed init!(s; prn=false)
        prn && @info "Loaded $name model in $time seconds"
        return s
    end
    sam = create_model("ram")
    tether_sam = create_model("tether")
    simple_sam = create_model("simple_ram")
    one_seg_sam = create_model("ram"; segments=1)
    one_seg_tether_sam = create_model("tether"; segments=1)
    return sam, tether_sam, simple_sam, one_seg_sam, one_seg_tether_sam
end

@setup_workload begin
    path = dirname(pathof(@__MODULE__))
    set_data_path(joinpath(path, "..", "data"))

    @compile_workload begin
        m1 = "Manifest-v1.$(VERSION.minor).toml"
        m2 = "Manifest-v1.$(VERSION.minor).toml.default"
        if filecmp(m1, m2)
            @info "Manifest files match, using the default xz files will work!"
            for input_path in readdir("data", join=true)
                if endswith(input_path, ".xz") && startswith("model", input_path)
                    output_path = replace(input_path, ".xz" => "")
                    decompress_binary(input_path, output_path)
                    println("Decompressed $input_path -> $output_path")
                end
            end
        else
            @warn "Manifest files differ, precompilation might be slow."
        end

        prn=true
        sam, tether_sam, simple_sam, _, _ = create_default_models(; prn)

        init!(sam; prn=false, reload=true)
        init!(sam; prn=false, reload=false)
        sim_oscillate!(sam; total_time=1.0)
        @show sam.sys_struct.name
        copy_to_simple!(sam, tether_sam, simple_sam; prn=false)
        find_steady_state!(sam)
        ss = SysState(sam)
        next_step!(sam)
        update_sys_state!(ss, sam)
        nothing
    end
end   
  
