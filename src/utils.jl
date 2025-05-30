function _propKepTE(r⃗::T, v⃗::T, Δt::Real, μ::Float64; itermax::Int=30, tol::Float64=1e-12)::Tuple{T, T} where T<:SVector

    # Defining Constants
    # @info "testing" x⃗₀ x⃗₀[1:3]
    # r⃗      = x⃗₀[1:3] |> SVector{3, eltype(x⃗₀)}; 
    r = norm(r⃗) 
    # v⃗      = x⃗₀[4:6] |> SVector{3, eltype(x⃗₀)};
    v = norm(v⃗)
    ε      = (0.5*v^2 - μ/r)
    a      = -0.5*μ/ε 
    ra, rμ = sqrt(a), sqrt(μ)
    σ      = dot(r⃗,v⃗)/rμ # Frequent Constant

    # Root Finding for Change in Eccentric Anomaly
    ΔE̅     =  Δt*sqrt(μ/a^3) # Vertical offset
    ΔE      = copy(ΔE̅)       # Initial guess
    for k in 1:itermax
        # Finding values for Halley method
        f   = -ΔE̅ + ΔE + (σ/ra)*(1 - cos(ΔE)) - (1 - r/a)*sin(ΔE)
        Δf  = 1 + (σ/ra)*sin(ΔE) + (r/a - 1)*cos(ΔE)
        Δ²f = (σ/ra)*cos(ΔE) - (r/a - 1)*sin(ΔE)

        # Solving and taking
        δE = -2*f*Δf/( 2*Δf^2 - f*Δ²f )
        ΔE += δE
        if abs(δE) < tol
            break
        end
    end

    # Creating Another Useful Constant
    cE, sE = cos(ΔE), sin(ΔE)
    ρ = a + (r - a)*cE + σ*(ra)*sE

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cE)
    g = a * σ / rμ * (1 - cE) + r * ra/rμ * sE
    ḟ = -rμ*ra / (ρ * r) * sE
    ġ = 1 - a / ρ * (1 - cE)

    # Converting to Cartesian States
    # r⃗ₖ = f*r⃗ + g*v⃗ # Most allocations
    # v⃗ₖ = ḟ*r⃗ + ġ*v⃗
    r⃗ₖ = @. f*r⃗ + g*v⃗ # Better
    v⃗ₖ = @. ḟ*r⃗ + ġ*v⃗
    
    return r⃗ₖ, v⃗ₖ
end

function _propKepTH(r⃗::T, v⃗::T, Δt::Real, μ::Float64; itermax::Int=30, tol::Float64=1e-12)::Tuple{T, T} where T<:SVector

    # Defining Constants
    r = norm(r⃗) 
    v = norm(v⃗)
    ε = (0.5*v^2 - μ/r)
    a = -0.5*μ/ε 
    σ = dot(r⃗,v⃗)/sqrt(μ) # Frequent Constant

    # Root Finding for Change in Hyperbolic Anomaly
    ΔH̅ = -Δt*sqrt(μ)/(-a)^(3/2) # Vertical offset
    ΔH = sign(Δt)               # Initial guess
    for k in 1:itermax
        # Finding values for Halley method
        f   = ΔH̅ - ΔH + (σ/sqrt(-a))*(cosh(ΔH) - 1) + (1 - r/a)*sinh(ΔH)
        Δf  = (σ/sqrt(-a))*sinh(ΔH) + (1 - r/a)*cosh(ΔH) - 1
        Δ²f = (σ/sqrt(-a))*cosh(ΔH) + (1 - r/a)*sinh(ΔH)

        # Solving and taking
        δH = -2*f*Δf/( 2*Δf^2 - f*Δ²f )
        ΔH += δH
        if abs(δH) < tol
            break
        end
    end
    # while abs(δH) > 1e-10
    #     δH = -func(ΔH)/∇func(ΔH)
    #     ΔH += δH
    # end


    # Creating Another Useful Constant
    cH, sH = cosh(ΔH), sinh(ΔH)
    ρ = a + (r - a)*cH + σ*sqrt(-a)*sH

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cH)
    g = a * σ / sqrt(μ) * (1 - cH) + r * sqrt(-a)/sqrt(μ) * sH
    ḟ = -sqrt(-a)*sqrt(μ) / (ρ * r) * sH
    ġ = 1 - a / ρ * (1 - cH)

    # Converting to Cartesian States
    r⃗ₖ = @. f*r⃗ + g*v⃗
    v⃗ₖ = @. ḟ*r⃗ + ġ*v⃗
    
    return r⃗ₖ, v⃗ₖ
end