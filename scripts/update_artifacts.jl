using ArtifactUtils
using HTTP
using JSON

# Get information about current release from GitHub API
function get_release_assets()
    # Extract info from environment
    repo = ENV["GITHUB_REPOSITORY"]
    tag = ENV["GITHUB_REF_NAME"]
    token = ENV["GITHUB_TOKEN"]

    api_url = "https://api.github.com/repos/$repo/releases/tags/$tag"

    headers = Dict("Authorization" => "token $token", "User-Agent" => "ArtifactUtils Script")

    response = HTTP.get(api_url, headers)
    if response.status != 200
        error("Failed to fetch release info from GitHub API: ", String(response.body))
    end
    release = JSON.parse(String(response.body))

    # Return list of asset info dicts
    release["assets"]
end

# Base URL for assets downloads
function asset_download_url(asset)
    # GitHub API Assets URL
    asset["browser_download_url"]
end

function main()
    println("Updating Artifacts.toml with release assets...")

    assets = get_release_assets()

    artifact_file = "Artifacts.toml"

    # For each .xz asset, add or update artifact entry
    for asset in assets
        name = asset["name"]
        if endswith(name, ".xz")
            url = asset_download_url(asset)
            local_path = joinpath("data", name)  # Dummy local path for hashing - should have .xz file locally or skip hash (can redownload)

            # Since local files don't exist in CI after build, fetch checksum from GitHub?
            # Alternative approach: calculate hash locally during create_xz_files step and save a checksum file
            # For simplicity, here we will download the file to temp and calculate hash:

            println("Downloading $name for checksum...")
            tmpfile = tempname()
            HTTP.download(url, tmpfile)

            # Compute SHA256 hash
            hash = ArtifactUtils.sha256(tmpfile)

            println("Adding artifact '$name' with checksum $hash and URL $url...")
            add_artifact!(artifact_file, replace(name, "." => "_"), url, sha256=hash; force=true)

            rm(tmpfile)
        end
    end

    println("Artifacts.toml updated.")
end

main()
