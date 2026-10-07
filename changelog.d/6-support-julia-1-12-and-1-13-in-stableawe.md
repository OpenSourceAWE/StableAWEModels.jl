### Changed
- Supported Julia versions are now 1.12 and 1.13; Julia 1.11 support is dropped.
- `bin/install -y` installs without asking anything, using the active Julia if it is supported and otherwise the newest installed supported one. The supported versions are those with a tracked `Manifest-vX.Y.toml.default`.
- `bin/install` no longer adds a `jl` alias to `~/.bashrc` or `~/.zshrc`.
