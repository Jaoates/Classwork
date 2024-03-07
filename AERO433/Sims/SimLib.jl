
module PowerFunction
    using SatelliteToolbox
    using ReferenceFrameRotations
    using LinearAlgebra

    export SimParams,rk,calcPower,getMode

    const re = EARTH_EQUATORIAL_RADIUS/1000

    struct SolarParams
        # this structure captures parameters for the model. it is expected to grow
        solarArrArea::Float64
        solarArrEff::Float64
        solarConst::Float64
    end

    struct PayloadParams
        sarDraw::Vector{Float64} # average wattage during opperation
        hyperDraw::Vector{Float64}
    end

    struct ElectricalParams
        comsDraw::Vector{Float64}
        obcDraw::Float64
    end

    struct GncParams
        gncDraw::Float64
    end

    struct PropsParams
        heaterDraw::Vector{Float64}
        engineDraw::Float64
    end

    struct SimParams
        solarparams::SolarParams
        payloadparams::PayloadParams
        electricalparams::ElectricalParams
        gncparams::GncParams
        propsparams::PropsParams
    end

    SolarParams() = SolarParams(3.5,.29,1320)
    PayloadParams() = PayloadParams([10,140],[10,600])
    ElectricalParams() = ElectricalParams([10,60],35)
    GncParams()=GncParams(42)
    PropsParams()=PropsParams([2,300],200)
    SimParams() = SimParams(SolarParams(),PayloadParams(),ElectricalParams(),GncParams(),PropsParams())


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

    function calcInShadow(r::AbstractVector,rsun::AbstractVector;radius::Number = re)
        # r is the spacecraft location in ECI
        # rsun is the location of the sun in ECI
        # radius is the radius of the body were orbiting
        # computes a simple calc based on the assumption that the sun is an infinitely far point source to tell if the s/c is in shadow.
        # assumes the shadow is a binary property

        # returns true if the spacecraft is in shadow /  eclipse

        if norm(r)<=radius
            throw(DomainError(r,"r is too small, norm(r)<radius"))
        end
        θA = 90 |> deg2rad
        θB = acos(radius/norm(r))
        θ = acos(dot(rsun,r)/(norm(rsun)*norm(r)))
        F = θA+θB<=θ # F = true if the sc is in shadow
    end

    function vecRotate(v,R::DCM)
        R*v
    end
    function vecRotate(v,R::Matrix)
        R*v
    end
    function vecRotate(v,R::Quaternion)
        v  = vect(R \ v * R)
    end

    function calcOverAfrica(jd::Number,X::AbstractVector)
        rECI = X[1:3]
        q = r_eci_to_ecef(Quaternion,J2000(),PEF(),jd)
        rECEF = vecRotate(rECI,q)
        Φ,λ,h=ecef_to_geodetic(rECEF*1000)
        Φbound = (-35.3,37.3).|> deg2rad # latitude
        λbound = (-17.2,52.6).|> deg2rad # longitude
        return Φbound[1]<Φ<Φbound[2] && λbound[1]<λ<λbound[2]
    end

    function calcSolarArrNormal(X::AbstractVector)
        r = X[1:3]
        n = r/norm(r)
        return n
    end

    function getMode(jd::Number,X::AbstractVector)
        r = X[1:3]

        mode=BitVector(undef,2)
        mode[1]=calcOverAfrica(jd::Number,X::AbstractVector)

        sun = sun_position_mod(jd)./1000
        Q = r_eci_to_eci(Quaternion,MOD(),J2000(),jd)
        sun = vecRotate(sun,Q)

        mode[2]=calcInShadow(r,sun)
        return mode
    end

    #
        
        # function Pcom(jd::Number,X::AbstractVector,params::SimParams)
        #     if calcOverAfrica(jd,X) return params.comsDraw else return 0 end
        # end

        # function Pobc(jd::Number,X::AbstractVector,params::SimParams)
        #     return params.obcDraw
        # end

        # function Pheater(jd::Number,X::AbstractVector,params::SimParams)
        #     r = X[1:3]
        #     sun = sun_position_mod(jd)./1000
        #     Q = r_eci_to_eci(Quaternion,MOD(),J2000(),jd)
        #     sun = vecRotate(sun,Q)

        #     if calcInShadow(r,sun) return params.heaterDraw else return 0 end
        # end
    #

    function Psun(jd::Number,X::AbstractVector,params::SimParams)
        r = X[1:3]
        v = X[4:6]
        A = params.solarparams.solarArrArea
        η = params.solarparams.solarArrEff
        Ks = params.solarparams.solarConst

        sun = sun_position_mod(jd)./1000
        Q = r_eci_to_eci(Quaternion,MOD(),J2000(),jd)
        sun = vecRotate(sun,Q)

        n = calcSolarArrNormal(X)
        p = dot(n,sun/norm(sun))*A*η*Ks

        if p<0 || calcInShadow(r,sun)
            p = 0
        end
        
        return p
    end

    function Psar(jd::Number,X::AbstractVector,params::SimParams)
        mode = getMode(jd,X)
        if mode[1] && mode[2] return params.payloadparams.sarDraw[2] else return params.payloadparams.sarDraw[1] end
    end

    function Phyper(jd::Number,X::AbstractVector,params::SimParams)
        mode = getMode(jd,X)
        if mode[1] && !mode[2] return params.payloadparams.hyperDraw[2] else return params.payloadparams.hyperDraw[1] end
    end

    function Pheater(jd::Number,X::AbstractVector,params::SimParams)
        mode = getMode(jd,X)
        if mode[2] return params.propsparams.heaterDraw[2] else return params.propsparams.heaterDraw[1] end
    end

    function Pcoms(jd::Number,X::AbstractVector,params::SimParams)
        mode = getMode(jd,X)
        if mode[1] return params.electricalparams.comsDraw[2] else return params.propsparams.heaterDraw[1] end
    end

    function Peng()
        # evertything payload is off
        # as much as 48 min per day
        # after GNC rotation
        # 
    end

    function calcPower(jd::Number,X::AbstractVector,params::SimParams)
        mode = getMode(jd,X)
        p = [
        Psun(jd,X,params),
        -Psar(jd,X,params),
        -Phyper(jd,X,params),
        -Pheater(jd,X,params),
        -Pcoms(jd,X,params),
        -params.electricalparams.obcDraw,
        -params.gncparams.gncDraw]
        return p
    end


    # function rk(r::AbstractVector,v::AbstractVector,jd::Number,params::SimParams)
    #     X = vcat(r,v)
    #     return (jd,X,params)
    # end



end


