using LambertsProblem
using StaticArrays, LinearAlgebra
using SPICE, Downloads

# =====================================================================
# === Setup

# Inputs
origin     = "Earth Barycenter"
target     = "Mars Barycenter"
launchSpan = ["2026-03-01T12:00:00", "2027-07-01T12:00:00"]
# launchSpan = ["2005-04-30T00:00:00", "2005-10-07T00:00:00"]
tofSpan    = [90.0, 1600.0]
tstep      = 1.0 # days
frame      = "ECLIPJ2000"
center     = "Solar System Barycenter"
gm         = 1.327111e11

# Retrieving and loading NAIF Kernels
kernels_files = [
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc",
]
kernels = [Downloads.download(f) for f in kernels_files]
for f in kernels; furnsh(f); end

# =====================================================================
# === Collecting Data

# Creating input ranges
lspan = utc2et(launchSpan[1]):tstep*86400:utc2et(launchSpan[2])
tspan = tofSpan[1]:tstep:tofSpan[2]

# Preallocating arrays
width, height = length(lspan), length(tspan)
v∞₊, v∞₋, sep, δν, R = zeros(width, height), zeros(width, height), zeros(width, height), zeros(width, height), zeros(width, height)

# Caculating Lamberts
for (k, launch) in enumerate(lspan)

    # Getting departure state
    dep = spkezr(origin, launch, frame, "None", center)[1]
    rₖ  = SA[dep[1:3]...]
    vₖ  = SA[dep[4:6]...]

    for (j, tof) in enumerate(tspan)

        # Getting arrival state
        arr = spkezr(target, launch + tof*86400, frame, "None", center)[1]
        rⱼ  = SA[arr[1:3]...]
        vⱼ  = SA[arr[4:6]...]

        # Creating lambert problem
        revs = 0
        prob = BallisticLambertsProblem(rₖ, rⱼ, tof*86400, gm, revs=0)
        sol  = solve(prob, Izzo(), silent=false)

        for rev_ in [1, 2]
            # Checking if multi-rev solution is better
            prob′ = BallisticLambertsProblem(rₖ, rⱼ, tof*86400, gm, revs=rev_)
            sol′ = solve(prob′, Izzo(), silent=true)
            if isvalid(sol′) && norm(sol.v⃗₁ - vₖ) > norm(sol′.v⃗₁ - vₖ)
                sol = sol′
                revs = rev_
                # @info "Hi" sol′
            end

            # Checking if multi-rev solution (longway) is better
            prob′ = BallisticLambertsProblem(rₖ, rⱼ, tof*86400, gm, revs=rev_, longway=true)
            sol′ = solve(prob′, Izzo(), silent=true)
            if isvalid(sol′) && norm(sol.v⃗₁ - vₖ) > norm(sol′.v⃗₁ - vₖ)
                sol = sol′
                revs = rev_
                # @info "Hi" sol′
            end
        end

        # Saving data
        v∞₊[k, j] = sol.v⃗₁ - vₖ |> norm
        v∞₋[k, j] = sol.v⃗₂ - vⱼ |> norm
        R[k, j]  = revs
        δν[k, j] = acosd( rₖ⋅rⱼ / norm(rₖ) / norm(rⱼ) )
        rₖ′ = SA[spkezr(origin, launch + tof*86400, frame, "None", center)[1][1:3]...]
        sep[k, j] = acosd( rₖ′⋅rⱼ / norm(rₖ′) / norm(rⱼ) )

    end

end

# =====================================================================
# === Plotting

using GLMakie, LaTeXStrings

glopts = (focus_on_show=true,)
GLMakie.activate!(inline=false; glopts...) # Display outside vscode
# GLMakie.activate!(inline=true; glopts...) # Display within vscode

begin
    fig = Figure()
    ax  = Axis(fig[1, 1], 
        title="Porkchop Plot",
        subtitle="$origin to $target",
        xlabel="Departure Date",
        # xticks=LinRange([utc2et(v) for v in launchSpan]..., 6),
        # Rounding date to nearest day, accounting for ET->UTC difference (69.1 seconds)
        xtickformat=vals->[replace(et2utc(round(v/86400)*86400+69.1, "ISOC", -3), "T"=>"\n") for v in vals],
        ylabel="Time of Flight [days]",
    )
    legend_handles = []
    legend_labels  = []

    # Adding revs used
    # c5 = contourf!(ax, lspan, tspan, R, 
    #     colormap=:Blues,
    #     # color=cgrad(:tab10).colors[6],
    #     # labels = true,
    #     # levels=0:30:180, 
    #     # depth_shift=-1000
    #     # transparency=0.1,
    #     inspectable=false,
    # )
    # push!(legend_handles, c5)
    # push!(legend_labels,  L"\textrm{Lambert Revs}")
    # Colorbar(fig[1, 2], c5, 
    #     ticks=[0, 1, 2],
    #     label=L"\textrm{Solution Revs}",
    # )

    # C3 Departure contours
    # c = contourf!(ax, lspan, tspan, v∞₊.^2, 
    #     # colormap=:coffee,
    #     colormap=:vik,
    #     # colormap=:seaborn_rocket_gradient,
    #     levels=LinRange(0.0, 50, 30)
    # )
    c = heatmap!(ax, lspan, tspan, v∞₊.^2, 
        # colormap=:coffee,
        colormap=:vik,
        # colormap=:seaborn_rocket_gradient,
        interpolate=true,
        colorrange=[0.0, 81.], 
        highclip=(:white, 0.0), # Making high end transparent
        inspector_label = (self, i, p) -> begin
            x, y = round(i[1]) |> Int, round(i[2]) |> Int
            launch = et2utc( lspan[x], "C", -3 )
            arrival = et2utc( lspan[x] + tspan[y]*86400, "C", -3  )
            """
            Launch: $launch
            Arrival: $arrival
            C₃: $(round(v∞₊[x, y]^2, digits=2)) km²/sec²
            Arrival v∞: $(round(v∞₋[x, y], digits=2)) km/sec
            Revs: $(R[x, y])
            """
        end,
    )
    contour!(ax, lspan, tspan, v∞₊.^2, 
        color=:black,
        labels=true,
        levels=0:5:80, 
        labelsize=12,
        # labelfont=:bold,
        inspectable=false,
    )
    Colorbar(fig[1, 2], c, 
        ticks=collect(3:12).^2,
        label=L"Launch $C_3$ $\left( km^2/s^2 \right)$",
    )

    # Adding vinf arrival contours
    c2 = contour!(ax, lspan, tspan, v∞₋, 
        color=cgrad(:tab10).colors[7],
        labels=true,
        levels=0:8, 
        labelsize=12,
        # labelfont=:bold,
        inspectable=false,
    )
    push!(legend_handles, c2)
    push!(legend_labels,  L"Arrival $v_\infty$")

    # Adding sep angle at arrival
    # c3 = contour!(ax, lspan, tspan, sep, 
    #     color=cgrad(:tab10).colors[5],
    #     labels = true,
    #     levels=0:30:180, 
        # inspectable=false,
    # )
    # push!(legend_handles, c3)
    # push!(legend_labels,  L"\textrm{Arrival SEP Angle}")

    # Adding transfer angle
    # c4 = contour!(ax, lspan, tspan, δν, 
    #     color=cgrad(:tab10).colors[3],
    #     labels = true,
    #     levels=0:30:180, 
        # inspectable=false,
    # )
    # push!(legend_handles, c4)
    # push!(legend_labels,  L"\textrm{Transfer Angle}")

    # Setting axis decorations
    tightlimits!(ax)
    # axislegend(ax, legend_handles, legend_labels,
    #     position=:lt,
    # )

    # Adding data inspector
    DataInspector(fig)

    # New event
    register_interaction!(ax, :click_index) do ev::MouseEvent, axis
        # Only react on the discrete “leftclick” event
        if ev.type === MouseEventTypes.leftclick
            # Reading event location
            launch, tof = ev.data
            launch_str = et2utc( launch, "C", -3 )
            arrival_str = et2utc( launch + tof*86400, "C", -3  )
            @info "Testing" launch_str arrival_str

            # Getting index from data position
            x = clamp(searchsortedfirst(lspan, launch),  1, length(lspan))
            y = clamp(searchsortedfirst(tspan, tof),  1, length(tspan))
            # @info "Testing" x y

            # Calculating values for lambert
            rₖ  = SA[spkezr(origin, launch,  frame, "None", center)[1][1:3]...]
            rⱼ  = SA[spkezr(target, launch + tof*86400, frame, "None", center)[1][1:3]...]
            # @info "Testing 2" rₖ rⱼ

            # # Solving lambert
            prob = BallisticLambertsProblem(rₖ, rⱼ, tof*86400, gm, revs=Int(R[x, y]))
            sol  = solve(prob, Izzo(), silent=false)
            # @info "Testing 3" sol

            # Plotting solution in inset plot    
            if isdefined(Main, :ax_inset); delete!(ax_inset); end
            ax_inset = Axis(fig[1, 1],
                backgroundcolor=:white,
                width=Relative(0.3),
                height=Relative(0.3),
                halign=0.1,
                valign=0.1,
                title="Lambert Arc",
                visible=false,
            )
            lines!(ax_inset, 1:500, cumsum(randn(500) .+ 0.5))

        end
    end

    # Setting up event (kinda works)
    #=
    on(events(ax).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            # Reading event location
            launch, tof = mouseposition(ax)
            launch_str = et2utc( launch, "C", -3 )
            arrival_str = et2utc( launch + tof*86400, "C", -3  )
            @info "Testing" launch_str arrival_str

            # Getting index from data position
            x = clamp(searchsortedfirst(lspan, launch),  1, length(lspan))
            y = clamp(searchsortedfirst(tspan, tof),  1, length(tspan))
            # @info "Testing" x y

            # Calculating values for lambert
            rₖ  = SA[spkezr(origin, launch,  frame, "None", center)[1][1:3]...]
            rⱼ  = SA[spkezr(target, launch + tof*86400, frame, "None", center)[1][1:3]...]
            # @info "Testing 2" rₖ rⱼ

            # # Solving lambert
            prob = BallisticLambertsProblem(rₖ, rⱼ, tof*86400, gm, revs=Int(R[x, y]))
            sol  = solve(prob, Izzo(), silent=false)
            # @info "Testing 3" sol

            # Plotting solution in inset plot
            ax_inset = Axis(fig[1, 1],
                backgroundcolor=:white,
                width=Relative(0.3),
                height=Relative(0.3),
                halign=0.1,
                valign=0.1,
                title="Lambert Arc",
            )
            lines!(ax_inset, 1:500, cumsum(randn(500) .+ 0.5))
                
        end
        return Consume(false)
    end
    =#

    fig
end
