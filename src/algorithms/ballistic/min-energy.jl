# =====================================================================
# === Hooking into interfaces

function solve(prob::MinimumEnergyLambertsProblem{T}) where T<:Real

    # Getting solution
    Δt,  _, v⃗₁, _, v⃗₂ = _min_energy(prob.r⃗₁, prob.r⃗₂, prob.μ, prob.revs, prob.retrograde)

    # Passing to output structure
    # (r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, Δtₑ, μ, revs, retrograde, longway)
    return MinimumEnergyLambertsSolution(prob.r⃗₁, v⃗₁, prob.r⃗₂, v⃗₂, Δt, prob.μ, revs=prob.revs, retrograde=prob.retrograde)
end


# =====================================================================
# === Solver

function _min_energy(
    r⃗1::SVector{3, <:Real}, 
    r⃗2::SVector{3, <:Real}, 
    μ::Real,
    revs::Int,
    retrograde::Bool,
    # longway::Bool,
    )

    # TODO: Catch for hyperbolic case
    # TODO: I don't think retrograde works either...

    # Processing inputs
    r1 = norm(r⃗1)
    r2 = norm(r⃗2)

    # Finding geometric parameters
    c⃗  = r⃗2 - r⃗1; c = norm(c⃗)                       # Chord
    δν = (1 - 2*retrograde)*acos( r⃗1⋅r⃗2 / r1 / r2 ) # Angle travelled
    s  = 0.5*(r1 + r2 + c)                          # Semi-perimeter

    # Solving for minimum energy ellipse
    aₘ = s/2                       # Semi-major axis
    pₘ = ((r1*r2)/c)*(1 - cos(δν)) # Semi-latus rectum
    # eₘ = sqrt( 1 - 2*pₘ/s )        # Eccentricity
    # α  = 2*asin(sqrt( s/2/aₘ ))    # Will always be π for minimum energy
    β  = 2*asin(sqrt( (s - c)/s ))

    # Solving for parameters of transfer
    Δt = sqrt(aₘ^3 / μ)*( 2π*revs + π + (β - sin(β)))
    # Δt = sqrt(aₘ^3 / μ)*( 2π*revs + π + (1 - 2*longway)*(β - sin(β)))
    f = 1 - r2/pₘ*(1 - cos(δν))
    g = r1*r2*sin(δν)/sqrt(μ*pₘ)
    ḟ = sqrt(μ/pₘ)*( (1 - cos(δν))/sin(δν) )*( (1 - cos(δν))/pₘ - 1/r1 - 1/r2 )
    ġ = 1 - r1/pₘ*(1 - cos(δν))

    # Creating outputs
    v⃗1 = (r⃗2 - f*r⃗1)/g
    v⃗2 = ḟ*r⃗1 + ġ*v⃗1

    return Δt, r⃗1, v⃗1, r⃗2, v⃗2

end



#=
if nargin < 5; shortlong = -1; end
scalarFlag = length(r1vec) == 1;

if ~scalarFlag
    r1 = norm(r1vec);
    r2 = norm(r2vec);
    r1dr2 = dot(r1vec, r2vec);
    dnu = atan2d(norm(cross(r1vec, r2vec)), dot(r1vec, r2v2c));
else
    r1 = r1vec;
    r2 = r2vec;
    r1dr2 = r1*r2*cosd(dnu);
end

% PAPER WAY
% pmin = ( (r1*r2)/sqrt(r1^2 + r2^2 - 2*r1dr2) )*(1 - cosd(dnu));
% 
% v1 = ( sqrt(mu*pmin)/(r1*r2*sind(dnu)) )*( r2vec - (1 - (r2/pmin)*(1 - cosd(dnu)))*r1vec );
% 
% amin = -0.5*mu / (0.5*norm(v1)^2 - mu/r1)
% 
% emin = sqrt(1 - pmin/amin)

% VALLADO
c = sqrt(r1^2 + r2^2 - 2*r1dr2);
s = 0.5*(r1 + r2 + c);
amin = s/2;
pmin = ((r1*r2)/c)*(1 - cosd(dnu));
emin = sqrt( 1 - 2*pmin/s );

beta = 2*asin(sqrt((s - c)/s));

tof = sqrt(amin^3 / mu)*(pi + shortlong*(beta - sin(beta)));

v1 = ( sqrt(mu*pmin)/(r1*r2*sind(dnu)) )*( r2vec - (1 - (r2/pmin)*(1 - cosd(dnu)))*r1vec );

lam = struct('v1', v1, 'amin', amin, 'pmin', pmin, 'emin', emin, 'tof', tof);
=#