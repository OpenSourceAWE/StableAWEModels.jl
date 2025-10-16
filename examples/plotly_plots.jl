# --------- imports ----------
using PlotlyJS          # direct PlotlyJS usage (no Plots.jl needed)
using PlotlyBase        # frame constructors & utilities
using SymbolicAWEModels # for STATIC/DYNAMIC tag
using LinearAlgebra

# --------- helpers ----------
# pick world pos if present, else CAD pos
get_pos(p) = hasproperty(p, :pos_w) ? p.pos_w : p.pos_cad



function plot3d_saddle(points, segments; title::AbstractString="3D Plotly Plot")
    # Use pos_cad if pos_w is not available
    get_pos(p) = hasproperty(p, :pos_w) ? p.pos_w : p.pos_cad
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]
    plt = Plots.scatter3d(x, y, z; markersize=2, markerstrokewidth=0, title=title, xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", legend=false)
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0
            @warn "Segment with zero index detected: $(s)"
            continue
        end
        Plots.plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]]; alpha=1, linewidth=1, color=:black)
    end
    fixed_idx = [i for (i, pt) in enumerate(points) if pt.type == SymbolicAWEModels.STATIC]
    !isempty(fixed_idx) && Plots.scatter3d!(x[fixed_idx], y[fixed_idx], z[fixed_idx]; markersize=2, markercolor=:red, markerstrokewidth=0)
    display(plt)
end



# --------- 1) one-off snapshot plot (your version) ----------
"""
    plot3d_v3(points, segments; title="3D Structure")

Create an interactive 3D snapshot of the structure with points and segments.
Fixed points in red, dynamic in blue, segments in black.
"""
function plot3d_v3(points, segments; title::AbstractString="3D Structure")
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]

    fixed_idx   = [i for (i, pt) in enumerate(points) if pt.type == SymbolicAWEModels.STATIC]
    dynamic_idx = [i for (i, pt) in enumerate(points) if pt.type != SymbolicAWEModels.STATIC]

    traces = PlotlyJS.GenericTrace[]

    # segments
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0
            @warn "Segment with zero index detected: $(s)"
            continue
        end
        push!(traces, PlotlyJS.scatter3d(
            x=[x[i], x[j]], y=[y[i], y[j]], z=[z[i], z[j]],
            mode="lines",
            line=PlotlyJS.attr(color="black", width=2),
            showlegend=false, hoverinfo="skip"
        ))
    end

    # dynamic points
    if !isempty(dynamic_idx)
        push!(traces, PlotlyJS.scatter3d(
            x=x[dynamic_idx], y=y[dynamic_idx], z=z[dynamic_idx],
            mode="markers",
            marker=PlotlyJS.attr(size=4, color="blue", symbol="circle"),
            name="Dynamic Points",
            hovertemplate="Point %{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>",
            text=string.(dynamic_idx)
        ))
    end

    # fixed points
    if !isempty(fixed_idx)
        push!(traces, PlotlyJS.scatter3d(
            x=x[fixed_idx], y=y[fixed_idx], z=z[fixed_idx],
            mode="markers",
            marker=PlotlyJS.attr(size=6, color="red", symbol="diamond"),
            name="Fixed Points",
            hovertemplate="Point %{text} (FIXED)<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>",
            text=string.(fixed_idx)
        ))
    end

    layout = PlotlyJS.Layout(
        title=title,
        scene=PlotlyJS.attr(
            xaxis=PlotlyJS.attr(title="X (m)"),
            yaxis=PlotlyJS.attr(title="Y (m)"),
            zaxis=PlotlyJS.attr(title="Z (m)"),
            aspectmode="data",
        ),
        showlegend=true,
        hovermode="closest"
    )

    plt = PlotlyJS.plot(traces, layout)
    display(plt)
    return plt
end

# --------- 2) animated version with play button and slider ----------
"""
    make_animated_plot3d(all_states, segments; title="…", dt=0.02)

Create an animated 3D plot with play button and slider.
`all_states` should be a vector of point position snapshots (each snapshot is a vector of points).
Returns the plot object.
"""
function make_animated_plot3d(all_states, segments; title::AbstractString="3D Structure", dt::Float64=0.02)
    n_frames = length(all_states)
    
    if n_frames == 0
        error("No states to animate")
    end
    
    # Build frames
    frames = PlotlyBase.PlotlyFrame[]
    frame_names = String[]
    
    for (frame_idx, points) in enumerate(all_states)
        x = [get_pos(p)[1] for p in points]
        y = [get_pos(p)[2] for p in points]
        z = [get_pos(p)[3] for p in points]
        
        # Build segment traces for this frame
        traces = PlotlyJS.GenericTrace[]
        for s in segments
            i, j = s.point_idxs
            if i == 0 || j == 0; continue; end
            push!(traces, PlotlyJS.scatter3d(
                x=[x[i], x[j]], y=[y[i], y[j]], z=[z[i], z[j]],
                mode="lines",
                line=PlotlyJS.attr(color="black", width=2),
                showlegend=false, hoverinfo="skip"
            ))
        end
        
        # Create PlotlyFrame - pass traces directly as first argument
        frame_name = "frame_$(frame_idx)"
        push!(frames, PlotlyBase.PlotlyFrame(PlotlyJS.attr(
            name = frame_name,
            data = traces,
            layout = PlotlyJS.attr(title = "$(title) - Step $(frame_idx)/$(n_frames)")
        )))
        push!(frame_names, frame_name)
    end
    
    # Initial trace (first frame)
    points = all_states[1]
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]
    
    initial_traces = PlotlyJS.GenericTrace[]
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0; continue; end
        push!(initial_traces, PlotlyJS.scatter3d(
            x=[x[i], x[j]], y=[y[i], y[j]], z=[z[i], z[j]],
            mode="lines",
            line=PlotlyJS.attr(color="black", width=2),
            showlegend=false, hoverinfo="skip"
        ))
    end
    
    # Layout with animation controls
    layout = PlotlyJS.Layout(
        title="$(title) - Step 1/$(n_frames)",
        scene=PlotlyJS.attr(
            xaxis=PlotlyJS.attr(title="X (m)"),
            yaxis=PlotlyJS.attr(title="Y (m)"),
            zaxis=PlotlyJS.attr(title="Z (m)"),
            aspectmode="data",
        ),
        showlegend=false,
        updatemenus=[
            PlotlyJS.attr(
                type="buttons",
                showactive=false,
                buttons=[
                    PlotlyJS.attr(
                        label="▶ Play",
                        method="animate",
                        args=[nothing, PlotlyJS.attr(
                            frame=PlotlyJS.attr(duration=dt*1000, redraw=true),
                            fromcurrent=true,
                            mode="immediate",
                            transition=PlotlyJS.attr(duration=0)
                        )]
                    ),
                    PlotlyJS.attr(
                        label="⏸ Pause",
                        method="animate",
                        args=[[nothing], PlotlyJS.attr(
                            frame=PlotlyJS.attr(duration=0, redraw=false),
                            mode="immediate",
                            transition=PlotlyJS.attr(duration=0)
                        )]
                    )
                ],
                x=0.1, y=0, xanchor="right", yanchor="top"
            )
        ],
        sliders=[
            PlotlyJS.attr(
                active=0,
                steps=[
                    PlotlyJS.attr(
                        args=[[frame_names[i]], PlotlyJS.attr(
                            frame=PlotlyJS.attr(duration=0, redraw=true),
                            mode="immediate",
                            transition=PlotlyJS.attr(duration=0)
                        )],
                        label="$(i)",
                        method="animate"
                    ) for (i, f) in enumerate(frames)
                ],
                x=0.1, y=0, len=0.9,
                xanchor="left", yanchor="top",
                pad=PlotlyJS.attr(b=10, t=50)
            )
        ]
    )
    
    # Create plot with frames as positional argument (third argument)
    plt = PlotlyJS.plot(initial_traces, layout, frames)
    display(plt)
    return plt
end

# --------- 3) simple live-updating version (no animation, just refresh) ----------
"""
    make_plot3d(points, segments; title="…") -> plt

Create a simple 3D plot showing only segments (no point markers).
Returns the plot object for manual updating if needed.
"""
function make_plot3d(points, segments; title::AbstractString="3D Structure")
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]

    traces = PlotlyJS.GenericTrace[]

    # Build segment traces
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0
            @warn "Segment with zero index detected: $(s)"
            continue
        end
        push!(traces, PlotlyJS.scatter3d(
            x=[x[i], x[j]], y=[y[i], y[j]], z=[z[i], z[j]],
            mode="lines",
            line=PlotlyJS.attr(color="black", width=2),
            showlegend=false, hoverinfo="skip"
        ))
    end

    layout = PlotlyJS.Layout(
        title=title,
        scene=PlotlyJS.attr(
            xaxis=PlotlyJS.attr(title="X (m)"),
            yaxis=PlotlyJS.attr(title="Y (m)"),
            zaxis=PlotlyJS.attr(title="Z (m)"),
            aspectmode="data",
        ),
        showlegend=false,
    )

    plt = PlotlyJS.plot(traces, layout)
    display(plt)
    return plt
end
