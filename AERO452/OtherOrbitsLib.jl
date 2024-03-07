# this place is to store things like curtis or Dr. A code
include("JoshOrbitsLib.jl")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function atmosphere(z)
    #
    # ATMOSPHERE calculates density for altitudes from sea level
    # through 1000 km using exponential interpolation.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##Geometric altitudes (km):
    h = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000]

    ##Corresponding densities (kg/m^3) from USSA76:
    r = [1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15]
    ##Scale heights (km):
    H = [ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 5.927 5.533 5.703 6.782 9.973 13.243 16.322 21.652 27.974 34.934 43.342 49.755 54.513 58.019 60.980 65.654 76.377 100.587 147.203 208.020]
    ##Handle altitudes outside of the range:
    if z > 1000
        z = 1000;
    elseif z < 0
        z = 0;
    end
    ##Determine the interpolation interval:
    for j = 1:27
        if z >= h[j] && z < h[j+1]
            i = j;
        end
    end

    if z == 1000
        i = 27;
    end
    ##Exponential interpolation:
    # println("curt")
    density = r[i]*exp(-(z - h[i])/H[i]);
end




# -------------------------------------------------------------------------
function lunar_position(jd)
#
#...Calculates the geocentric equatorial position vector of the moon
# given the Julian day.
#
# User M-functions required: None
# -------------------------------------------------------------------------
#...Earth's radius (km):
RE = re;
# ------------------------- implementation -----------------
#...Time in centuries since J2000:
T = (jd - 2451545)/36525;
#...Ecliptic longitude (deg):
e_long = 218.32 + 481267.881*T+ 6.29*sind(135.0 + 477198.87*T) - 1.27*sind(259.3 - 413335.36*T)+ 0.66*sind(235.7 + 890534.22*T) + 0.21*sind(269.9 + 954397.74*T)- 0.19*sind(357.5 + 35999.05*T) - 0.11*sind(186.5 + 966404.03*T);
e_long = mod(e_long,360);
#...Ecliptic latitude (deg):
e_lat = 5.13*sind( 93.3 + 483202.02*T) + 0.28*sind(228.2 + 960400.89*T)- 0.28*sind(318.3 + 6003.15*T) - 0.17*sind(217.6 - 407332.21*T);
e_lat = mod(e_lat,360);
#...Horizontal parallax (deg):
h_par = 0.9508+ 0.0518*cosd(135.0 + 477198.87*T) + 0.0095*cosd(259.3 - 413335.36*T)+ 0.0078*cosd(235.7 + 890534.22*T) + 0.0028*cosd(269.9 + 954397.74*T);
h_par = mod(h_par,360);
#...Angle between earth's orbit and its equator (deg):
obliquity = 23.439291 - 0.0130042*T;
#...Direction cosines of the moon's geocentric equatorial position vector:
l = cosd(e_lat) * cosd(e_long);
m = cosd(obliquity)*cosd(e_lat)*sind(e_long) - sind(obliquity)*sind(e_lat);
n = sind(obliquity)*cosd(e_lat)*sind(e_long) + cosd(obliquity)*sind(e_lat);
#...Earth-moon distance (km):
dist = RE/sind(h_par);
#...Moon's geocentric equatorial position vector (km):
r_moon = dist.*[l,m,n];
end #lunar_position
# -------------------------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function solar_position(jd)
#
# This function calculates the geocentric equatorial position vector
# of the sun, given the julian date.
#
# User M-functions required: None
# -------------------------------------------------------------------------
#...Astronomical unit (km):
AU = 149597870.691;
#...Julian days since J2000:
n = jd - 2451545;
#...Julian centuries since J2000:
cy = n/36525;
#...Mean anomaly (deg{:
M = 357.528 + 0.9856003*n;
M = mod(M,360);
#...Mean longitude (deg):
L = 280.460 + 0.98564736*n;
L = mod(L,360);
#...Apparent ecliptic longitude (deg):
lamda = L + 1.915*sind(M) + 0.020*sind(2*M);
lamda = mod(lamda,360);
#...Obliquity of the ecliptic (deg):
eps = 23.439 - 0.0000004*n;
#...Unit vector from earth to sun:
u = [cosd(lamda), sind(lamda)*cosd(eps), sind(lamda)*sind(eps)];
#...Distance from earth to sun (km):
rS = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*AU;
#...Geocentric position vector (km):
r_S = rS*u;
return (r_S,lamda,eps)
end #solar_position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# # test case for solar and lunar vectors
# jd = DateTime(1997,8,10)
# jd = date_to_jd(jd)
# norm(moon_position_mod(jd,Val(:Vallado))./1000-lunar_position(jd))
# norm(sun_position_mod(jd)./1000-solar_position(jd)[1])


# This is renos lamberts function
S(z) = sum(k -> ((-1)^k * z.^k./factorial(2*k+3)), 0:8)
C(z) = sum(k -> ((-1)^k * z.^k./factorial(2*k+2)), 0:8)
function lambert_reno(r₁::AbstractVector, r₂::AbstractVector, Δt::Number; pro = 1, tol = 1e-8, μ = μe)
    Δθ = acosd(dot(r₁, r₂)/(norm(r₁) * norm(r₂)))
    if pro * cross(r₁, r₂)[3] < 0
        Δθ = 360 - Δθ
    end
    
    A = sind(Δθ)*sqrt(norm(r₁)*norm(r₂) / (1 - cosd(Δθ)))
    @assert abs(A) > tol

    y(z) = norm(r₁) + norm(r₂) + A*(z*S(z) - 1)/sqrt(C(z))
    F(z) = (y(z)/C(z))^1.5 * S(z) + A * sqrt(y(z)) - sqrt(μ)*Δt
    z_bounds = (-4*pi^2, 4*pi^2)
    z = find_zero(F, 0, rtol=tol)

    f = 1 - y(z)/norm(r₁)
    g = A * sqrt(y(z)/μ)
    ḟ = sqrt(μ)/(norm(r₁)*norm(r₂)) * sqrt(y(z)/C(z)) * (z*S(z) - 1)
    ġ = 1 - y(z)/norm(r₂)

    v₁ = 1/g * (r₂ - f*r₁)
    v₂ = 1/g * (ġ*r₂ - r₁)
    return v₁, v₂
end