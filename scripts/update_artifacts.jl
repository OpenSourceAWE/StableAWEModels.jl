# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("JSON3" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using JSON3
using Downloads
using SHA
using Tar
using TOML
using CodecZlib

# Params
owner = "OpenSourceAWE"
repo = "SymbolicAWEModels"
project_toml_path = joinpath(@__DIR__, "..", "Project.toml")
artifacts_toml_path = joinpath(@__DIR__, "..", "Artifacts.toml")

# Extract version from Project.toml
function get_version(path)
    toml = TOML.parsefile(path)
    return toml["version"]
end

version = get_version(project_toml_path)
tag = "v$version"
println("Detected version: $version, using tag: $tag")

# Query GitHub API for release by tag
function get_release_by_tag(owner, repo, tag)
    url = "https://api.github.com/repos/$owner/$repo/releases/tags/$tag"
    resp = Downloads.download(url; headers=["Accept" => "application/vnd.github.v3+json"])
    open(resp) do io
        return JSON3.read(io)
    end
end

release = get_release_by_tag(owner, repo, tag)

# Download asset to temp file
function download_asset(url)
    temp_path = mktemp(suffix=".tar.gz")[1]
    Downloads.download(url, temp_path)
    return temp_path
end

# Compute sha256 hash of file
function sha256_of_file(path)
    open(path) do io
        return bytes2hex(sha256(io))
    end
end

# Compute git-tree-sha1 of decompressed tarball content
function git_tree_sha1_of_tarball(path)
    buf = IOBuffer()
    open(GzipDecompressorStream, path) do gz
        write(buf, read(gz))
    end
    seekstart(buf)
    return bytes2hex(Tar.tree_hash(buf))
end

# Process assets and generate Artifacts.toml content
artifacts = Dict{String,Any}()

for asset in release["assets"]
    fname = asset["name"]
    url = asset["browser_download_url"]
    if endswith(fname, ".tar.gz")
        println("Processing asset: $fname")
        local_path = download_asset(url)
        sha256hash = sha256_of_file(local_path)
        gtree_sha1 = git_tree_sha1_of_tarball(local_path)

        artifact_name = replace(fname, r"\.tar\.gz" => "")
        artifacts[artifact_name] = Dict(
            "git-tree-sha1" => gtree_sha1,
            "download" => [
                Dict(
                    "url" => url,
                    "sha256" => sha256hash
                )
            ]
        )
        rm(local_path; force=true)
    end
end

# Write Artifacts.toml
toml_dict = Dict(artifact => Dict(v) for (artifact,v) in artifacts)

open(artifacts_toml_path, "w") do io
    TOML.print(io, toml_dict)
end

println("Artifacts.toml successfully written to $artifacts_toml_path")
