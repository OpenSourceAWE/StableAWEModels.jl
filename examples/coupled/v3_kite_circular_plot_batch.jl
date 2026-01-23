using CSV
using DataFrames
using Statistics
using Printf
using GLMakie
using Colors: distinguishable_colors

function load_batch_csv(batch_name::AbstractString)
    batch_dir = joinpath("processed_data", "v3_kite", batch_name)
    csv_path = joinpath(batch_dir, "circle_batch_analysis.csv")
    isfile(csv_path) || error("CSV not found: $csv_path")
    df = CSV.read(csv_path, DataFrame)
    return df, batch_dir
end

function plot_batch(df::DataFrame; batch_dir::AbstractString)
    df = filter(row -> isfinite(row.us) && isfinite(row.v_app) &&
                       isfinite(row.yaw_rate) && isfinite(row.yaw_rate_paper) &&
                       isfinite(row.up), df)
    if :lt in names(df)
        df = filter(row -> isfinite(row.lt), df)
    end

    ups = sort(unique(df.up))
    palette = distinguishable_colors(length(ups))

    fig = Figure(size=(800, 350))
    ax1 = Axis(fig[1, 1], xlabel="us*v_app [m/s]", ylabel="yaw_rate_paper [deg/s]")
    ax2 = Axis(fig[1, 2], xlabel="us [-]", ylabel="CS [-]")

    legend_handles = Any[]
    legend_labels = String[]
    plotted_any = false

    for (i, up_val) in enumerate(ups)
        rows = df[df.up .== up_val, :]
        x = rows.us .* rows.v_app

        mask1 = isfinite.(x) .& isfinite.(rows.yaw_rate_paper)
        mask2 = isfinite.(rows.us) .& isfinite.(rows.cs)

        color = palette[i]
        label = @sprintf("up=%.3f", up_val)

        # Report mean percent delta between axis 1 and axis 2 data (shared points)
        mask12 = mask1 .& mask2
        if any(mask12)
            x1 = x[mask12]
            x2 = x[mask12]
            y1 = rows.yaw_rate[mask12]
            y2 = rows.yaw_rate_paper[mask12]

            nz_y = abs.(y2) .> 1e-12
            pct_y = (y1[nz_y] .- y2[nz_y]) ./ y2[nz_y] .* 100

            nz_x = abs.(x2) .> 1e-12
            pct_x = (x1[nz_x] .- x2[nz_x]) ./ x2[nz_x] .* 100

            @info "Mean percent delta (yaw_rate vs yaw_rate_paper)" up=up_val mean_pct=isempty(pct_y) ? NaN : mean(pct_y) n=length(pct_y)
            @info "Mean percent delta (x axis 1 vs 2)" up=up_val mean_pct=isempty(pct_x) ? NaN : mean(pct_x) n=length(pct_x)
        end

        if any(mask1)
            h = scatter!(ax1, x[mask1], rows.yaw_rate_paper[mask1];
                         color=color, markersize=8)
            push!(legend_handles, h)
            push!(legend_labels, label)
            plotted_any = true
        end

        if any(mask2)
            scatter!(ax2, rows.us[mask2], rows.cs[mask2];
                     color=color, markersize=8)
        end
    end

    # Plot dynamic time-series data points if available
    usva_cols = filter(x -> startswith(string(x), "usva_"), names(df))
    yaw_rate_cols = filter(x -> occursin(r"^yaw_rate_\d+$", string(x)), names(df))
    
    if !isempty(usva_cols) && !isempty(yaw_rate_cols)
        @info "Found time-series columns" n_usva=length(usva_cols) n_yaw_rate=length(yaw_rate_cols)
        
        # Plot dynamic data for each up value
        for (i, up_val) in enumerate(ups)
            rows = df[df.up .== up_val, :]
            color = palette[i]
            
            # Collect all dynamic points for this up value
            dynamic_x = Float64[]
            dynamic_y = Float64[]
            
            for row in eachrow(rows)
                for (usva_col, yaw_col) in zip(usva_cols, yaw_rate_cols)
                    x_val = row[usva_col]
                    y_val = row[yaw_col]
                    if isfinite(x_val) && isfinite(y_val)
                        push!(dynamic_x, x_val)
                        push!(dynamic_y, y_val)
                    end
                end
            end
            
            if !isempty(dynamic_x)
                label = i == 1 ? "dynamic" : ""
                h = scatter!(ax1, dynamic_x, dynamic_y;
                             color=color,
                             markersize=4, marker=:circle, alpha=0.5)
                if !isempty(label)
                    push!(legend_handles, h)
                    push!(legend_labels, label)
                end
                plotted_any = true
            end
        end
    end

    if plotted_any && !isempty(legend_handles)
        axislegend(ax1, legend_handles, legend_labels; position=:rb)
    else
        @warn "No valid points to plot; legend skipped"
    end

    lt_tag = ""
    if :lt in names(df)
        lt_vals = unique(df.lt)
        if length(lt_vals) == 1
            lt_tag = "_lt_$(Int(round(lt_vals[1])))"
        end
    end
    out_path = joinpath(batch_dir, "circle_batch_plot$(lt_tag).png")
    save(out_path, fig)
    @info "Saved plot" path=out_path
end

function main()
    batch_name = isempty(ARGS) ? "" : strip(ARGS[1])
    batch_name = "circular_2025_batch_2026_01_11_11_29_19"
    if isempty(batch_name)
        print("Enter batch folder name (e.g. batch_2026_01_07_10_04_38): ")
        batch_name = strip(readline())
    end
    isempty(batch_name) && error("Batch folder name is required.")

    df, batch_dir = load_batch_csv(batch_name)
    plot_batch(df; batch_dir=batch_dir)
end

main()
