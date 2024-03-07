include("JoshOrbitsLib1.jl")
using Statistics

function axesEqual!(radius = re+500)
    plot!(xlimit = [-radius,radius],ylimit = [-radius,radius],zlimit = [-radius,radius])
end

function uvec(r)
    # gets the unit vector of the vector r
    r̂ = r/norm(r)
end

function vecAngle(r1,r2)
    # calculate the angle between two vectors
    acos(dot(r1,r2)/(norm(r1)*norm(r2)))
end

function calcRaDec(r)
    # takes r and returns

    # r (ECI)

    # right Ascension (rad)
    # Declination (rad)
    l,m,n = r |> uvec
    δ = asin(n)
    if m>0
        α = acos(l/cos(δ))
    else
        α = 2pi-acos(l/cos(δ))
    end
    return α,δ
end

function calcLocHourAngle(θ,α)
    # Local Sideral Time (radians)
    # Right Ascension (radians)

    # Local Hour Angle (Radians)

    LHA = mod(θ - α,deg2rad(360))
end

function calcLst(λ,jd)
    # geodectic longitude (radians)
    # julian date (days)

    # lst in rad
    λ = λ|>rad2deg
    T_UT1 = (jd-2451545)./36525;
    ThetaGMST = 67310.54841 + (876600*3600 + 8640184.812866).*T_UT1 + .093104 .*(T_UT1.^2) - (6.2*10^-6).*(T_UT1.^3);
    ThetaGMST = mod((mod(ThetaGMST,86400*(ThetaGMST./abs(ThetaGMST)))/240),360);
    θ = ThetaGMST + λ;
    θ = θ|>deg2rad
end



function calcElAz(Φ::Number,λ::Number,α::Number,δ::Number,jd::Number)
    # function and subfunctions are verified 🙄
        # see matlab script

    # geodectic latitude (radians)
    # geodectic longitude (radians)
    # right Ascension (radians)
    # declination (radians)
    # julian date (days)

    # appearant elevation (radians)
    # appearant azimuth (radians)

    lst = calcLst(λ,jd)
    ha = calcLocHourAngle(lst,α) 
    el = asin(sin(Φ).*sin(δ)+cos(Φ).*cos(δ).*cos(ha));
    az = mod(atan(-sin(ha).*cos(δ)./cos(el),(sin(δ)-sin(el).*sin(Φ))./(cos(el).*cos(Φ))),deg2rad(360));
    el = el
    az = az
    return el,az
end

#  ↓functionally equivalent to RAZEL as a whole
function calcElAz(Φ::Number,λ::Number,ρ::AbstractVector,jd::Number)
        α,δ = calcRaDec(ρ)
        return calcElAz(Φ::Number,λ::Number,α::Number,δ::Number,jd::Number) 
end

# function rAZEL(ρ,Φ,λ,lst)
#     # its walmart RAZEL
#     R3(θ) = angle_to_dcm(θ|>deg2rad,:Z)
#     R2(θ) = angle_to_dcm(θ|>deg2rad,:Y)
#     R1(θ) = angle_to_dcm(θ|>deg2rad,:X)

#     rho = ρ
#     lat = Φ |> rad2deg
#     lon = λ |> rad2deg
#     lst = lst |> rad2deg
#     rhoecef = R3(lst-lon)*rho;
#     rhosez = R2(90-lat)*R3(lon)*rhoecef
#     el = asin(rhosez[3]/norm(rhosez)); 
#     el = mod(el,2*pi);
#     if el != pi/2
#        if atan(rhosez[1],rhosez[2]) < pi/2
           
#            beta = acos(-rhosez[1]/sqrt(rhosez[1]^2 + rhosez[2]^2));
#            beta = mod(beta,2*pi);
#        else
           
#            beta = acos(-rhosez[1]/sqrt(rhosez[1]^2 + rhosez[2]^2));
#            beta = 2*pi - mod(beta,2*pi);
#        end
#     else
        
#         beta = 0;
#     end
#     el = el*180/pi;
#     beta = beta*180/pi;
#     return el|>deg2rad,beta|>deg2rad
# end


function rhoElAzToNED(rho,el,az)
    # rho as a distance
    # el in rad
    # az in rad

    # rho in NED
    rhosez = Vector{Float64}(undef,3) 

    rhosez[1] = -rho*cos(el)*cos(az);
    rhosez[2] = rho*cos(el)*sin(az);
    rhosez[3] = rho*sin(el);
    rhoNED = rhosez.*[-1,1,-1];
end

function rhoElAzToECI(rho,el,az,Φ,λ,h,jd)
    # rho as a distance
    # el in rad
    # az in rad
    # Φ in rad
    # λ in rad
    # h in km
    # jd in days

    rhoNED = rhoElAzToNED(rho,el,az)
    rhoECEF = ned_to_ecef(rhoNED,Φ,λ,h*1000)
    q = r_ecef_to_eci.(Quaternion,PEF(),J2000(),jd)
    rhoECI = vecRotate(rhoECEF,q)
end

function doHgibbs(r1,r2,r3,jd1,jd2,jd3;μ = μe,tolangle= 0.01745329251994)
    # function from Vallado


    theta = 0.0;
    theta1= 0.0;
    magr1 = norm( r1 );
    magr2 = norm( r2 );
    magr3 = norm( r3 );
    for i= 1 : 3
        v2(i)= 0.0;
      end

    
    dt21= (jd2-jd1)*86400.0;
    dt31= (jd3-jd1)*86400.0;    #% differences in times in secs
    dt32= (jd3-jd2)*86400.0;

    p = cross( r2,r3 );
    pn = uvec( p );
    r1n = uvec( r1 );
    copa=  asin( dot( pn,r1n ) );
    if ( abs( dot(r1n,pn) ) > tolangle )
        # error= "not coplanar";
        @warn "not coplanar"
    end

   theta  = vecAngle( r1,r2 );
    theta1 = vecAngle( r2,r3 );
    if ( (theta > tolangle) | (theta1 > tolangle) )  
        # error= "   angl > 1ø";\
        @warn "angl > 1ø"
    end

    term1= -dt32*( 1.0/(dt21*dt31) + μ/(12.0*magr1*magr1*magr1) );
    term2= (dt32-dt21)*( 1.0/(dt21*dt32) + μ/(12.0*magr2*magr2*magr2) );
    term3=  dt21*( 1.0/(dt32*dt31) + μ/(12.0*magr3*magr3*magr3) );

    v2 =  term1*r1 + term2* r2 + term3* r3;

    return (v2,theta,theta1,copa)
end
function doHgibbs(rs,jds;μ = μe,tolangle= 0.01745329251994)
    r1,r2,r3 = rs
    jd1,jd2,jd3 = jds
    return doHgibbs(r1,r2,r3,jd1,jd2,jd3;μ = μ,tolangle= tolangle)
end

function doWeightedLLSdealio(more,times,Rsite)
    # yes youll need to prove this works rather rigorously
    rs = calcrFromAngles()
    vs = GibbsLambertsETC()
    # propagate them to one of the times (pass in kwargs for which propagator?) yes this is repetative
    # average these values to get Xnom  (basiclly our best initial guess)
    
    # while RMSprev-RMScur/RMSprev < tol (1e-3)
    # need a new calcAzEl and or RA DEC plug Xnom in and get out what we think the observations should have between
    # y = observed - calced (above)
    # Xmod[i] = Xnom + δ[i]
    # δ[i] = .001*Xnom[i]
    # calc AzEl/RaDec again
    # stick these into their own collunm of H (2x6) or (3x6)
    # if this is the second or higher iter, add these values to H, but be sure to propagate this value to t1
    x̂ = inv(H'*w*H)*H'*w*y # check if SVD needed for the inv\
    Xnom +=x̂

    return NaN
end


# function calcRMS(data)
    #     # W is wiehgts
    #     # N is number observation
    #     # n dim of observation
    #     return RMS = sqrt(transpose(y)*W*y/(n*N))
    # end


# function calcrFromAngles(angles,times,Rsite)
    #     # takes az el or ra dec and gives rECI (or ECEF, idk which is better)
    #     return NaN
    # end
        
# function GibbsLambertsETC(rs,ts)
    #     # need 3 of each r,t
    #     # do Gibbsy poopy or pick another
    #     # propagate them to one of the times? (pass in kwargs for which propagator?)
    #     # return v
    #     return NaN
    # end

# function calcLagrangeCoeffs(alsoidontknow)
    #     [[f,g],[ḟ,ġ]] * NaN
    # end

# function calcStateTransMat(yoNoSe)
    #     calcLagrangeCoeffs()
    #     𝚽 = [[f,g],[ḟ,ġ]]* NaN
    # end

# function calcGainMatrix()
    #     # see W5L2 slide 14
    #     # covariance (hat tho)
    #     # H...
    #     # weights
    #     K = P̂*transpose(H)*W + inv(P̂) # needs a bit more alebra to be calcuable
    # end


