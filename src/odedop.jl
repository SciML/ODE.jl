type DOP853
    a21::Float64
    a31::Float64
    a32::Float64
    a41::Float64
    a43::Float64
    a51::Float64
    a53::Float64
    a54::Float64
    a61::Float64
    a64::Float64
    a65::Float64
    a71::Float64
    a74::Float64
    a75::Float64
    a76::Float64
    a81::Float64
    a84::Float64
    a85::Float64
    a86::Float64
    a87::Float64
    a91::Float64
    a94::Float64
    a95::Float64
    a96::Float64
    a97::Float64
    a98::Float64
    a101::Float64
    a104::Float64
    a105::Float64
    a106::Float64
    a107::Float64
    a108::Float64
    a109::Float64
    a111::Float64
    a114::Float64
    a115::Float64
    a116::Float64
    a117::Float64
    a118::Float64
    a119::Float64
    a1110::Float64
    a121::Float64
    a124::Float64
    a125::Float64
    a126::Float64
    a127::Float64
    a128::Float64
    a129::Float64
    a1210::Float64
    a1211::Float64
    b1::Float64
    b6::Float64
    b7::Float64
    b8::Float64
    b9::Float64
    b10::Float64
    b11::Float64
    b12::Float64
    c2::Float64
    c3::Float64
    c4::Float64
    c5::Float64
    c6::Float64
    c7::Float64
    c8::Float64
    c9::Float64
    c10::Float64
    c11::Float64
    bhh1::Float64
    bhh2::Float64
    bhh3::Float64
    er1::Float64
    er6::Float64
    er7::Float64
    er8::Float64
    er9::Float64
    er10::Float64
    er11::Float64
    er12::Float64
    c14::Float64
    c15::Float64
    c16::Float64
    d41::Float64
    d46::Float64
    d47::Float64
    d48::Float64
    d49::Float64
    d410::Float64
    d411::Float64
    d412::Float64
    d413::Float64
    d414::Float64
    d415::Float64
    d416::Float64
    d51::Float64
    d56::Float64
    d57::Float64
    d58::Float64
    d59::Float64
    d510::Float64
    d511::Float64
    d512::Float64
    d513::Float64
    d514::Float64
    d515::Float64
    d516::Float64
    d61::Float64
    d66::Float64
    d67::Float64
    d68::Float64
    d69::Float64
    d610::Float64
    d611::Float64
    d612::Float64
    d613::Float64
    d614::Float64
    d615::Float64
    d616::Float64
    d71::Float64
    d76::Float64
    d77::Float64
    d78::Float64
    d79::Float64
    d710::Float64
    d711::Float64
    d712::Float64
    d713::Float64
    d714::Float64
    d715::Float64
    d716::Float64
    a141::Float64
    a147::Float64
    a148::Float64
    a149::Float64
    a1410::Float64
    a1411::Float64
    a1412::Float64
    a1413::Float64
    a151::Float64
    a156::Float64
    a157::Float64
    a158::Float64
    a1511::Float64
    a1512::Float64
    a1513::Float64
    a1514::Float64
    a161::Float64
    a166::Float64
    a167::Float64
    a168::Float64
    a169::Float64
    a1613::Float64
    a1614::Float64
    a1615::Float64
end

type DOPRI5
    c2::Float64
    c3::Float64
    c4::Float64
    c5::Float64
    a21::Float64
    a31::Float64
    a32::Float64
    a41::Float64
    a42::Float64
    a43::Float64
    a51::Float64
    a52::Float64
    a53::Float64
    a54::Float64
    a61::Float64
    a62::Float64
    a63::Float64
    a64::Float64
    a65::Float64
    a71::Float64
    a73::Float64
    a74::Float64
    a75::Float64
    a76::Float64
    e1::Float64
    e3::Float64
    e4::Float64
    e5::Float64
    e6::Float64
    e7::Float64
    d1::Float64
    d3::Float64
    d4::Float64
    d5::Float64
    d6::Float64
    d7::Float64
end

const dopri5coeff = DOPRI5(
    0.2,
    0.3,
    0.8,
    8.0/9.0,
    0.2,
    3.0/40.0,
    9.0/40.0,
    44.0/45.0,
    -56.0/15.0,
    32.0/9.0,
    19372.0/6561.0,
    -25360.0/2187.0,
    64448.0/6561.0,
    -212.0/729.0,
    9017.0/3168.0,
    -355.0/33.0,
    46732.0/5247.0,
    49.0/176.0,
    -5103.0/18656.0,
    35.0/384.0,
    500.0/1113.0,
    125.0/192.0,
    -2187.0/6784.0,
    11.0/84.0,
    71.0/57600.0,
    -71.0/16695.0,
    71.0/1920.0,
    -17253.0/339200.0,
    22.0/525.0,
    -1.0/40.0  ,
    -12715105075.0/11282082432.0,
    87487479700.0/32700410799.0,
    -10690763975.0/1880347072.0,
    701980252875.0/199316789632.0,
    -1453857185.0/822651844.0,
    69997945.0/29380423.0)

const dop853coeff = DOP853(5.26001519587677318785587544488e-2,
    1.97250569845378994544595329183e-2,
    5.91751709536136983633785987549e-2,
    2.95875854768068491816892993775e-2,
    8.87627564304205475450678981324e-2,
    2.41365134159266685502369798665e-1,
    -8.84549479328286085344864962717e-1,
    9.24834003261792003115737966543e-1,
    3.7037037037037037037037037037e-2,
    1.70828608729473871279604482173e-1,
    1.25467687566822425016691814123e-1,
    3.7109375e-2,
    1.70252211019544039314978060272e-1,
    6.02165389804559606850219397283e-2,
    -1.7578125e-2,
    3.70920001185047927108779319836e-2,
    1.70383925712239993810214054705e-1,
    1.07262030446373284651809199168e-1,
    -1.53194377486244017527936158236e-2,
    8.27378916381402288758473766002e-3,
    6.24110958716075717114429577812e-1,
    -3.36089262944694129406857109825e0,
    -8.68219346841726006818189891453e-1,
    2.75920996994467083049415600797e1,
    2.01540675504778934086186788979e1,
    -4.34898841810699588477366255144e1,
    4.77662536438264365890433908527e-1,
    -2.48811461997166764192642586468e0,
    -5.90290826836842996371446475743e-1,
    2.12300514481811942347288949897e1,
    1.52792336328824235832596922938e1,
    -3.32882109689848629194453265587e1,
    -2.03312017085086261358222928593e-2,
    -9.3714243008598732571704021658e-1,
    5.18637242884406370830023853209e0,
    1.09143734899672957818500254654e0,
    -8.14978701074692612513997267357e0,
    -1.85200656599969598641566180701e1,
    2.27394870993505042818970056734e1,
    2.49360555267965238987089396762e0,
    -3.0467644718982195003823669022e0,
    2.27331014751653820792359768449e0,
    -1.05344954667372501984066689879e1,
    -2.00087205822486249909675718444e0,
    -1.79589318631187989172765950534e1,
    2.79488845294199600508499808837e1,
    -2.85899827713502369474065508674e0,
    -8.87285693353062954433549289258e0,
    1.23605671757943030647266201528e1,
    6.43392746015763530355970484046e-1,
    5.42937341165687622380535766363e-2,
    4.45031289275240888144113950566e0,
    1.89151789931450038304281599044e0,
    -5.8012039600105847814672114227e0,
    3.1116436695781989440891606237e-1,
    -1.52160949662516078556178806805e-1,
    2.01365400804030348374776537501e-1,
    4.47106157277725905176885569043e-2,
    0.526001519587677318785587544488e-01,
    0.789002279381515978178381316732e-01,
    0.118350341907227396726757197510e+00,
    0.281649658092772603273242802490e+00,
    0.333333333333333333333333333333e+00,
    0.25e+00,
    0.307692307692307692307692307692e+00,
    0.651282051282051282051282051282e+00,
    0.6e+00,
    0.857142857142857142857142857142e+00,
    0.244094488188976377952755905512e+00,
    0.733846688281611857341361741547e+00,
    0.220588235294117647058823529412e-01,
    0.1312004499419488073250102996e-01,
    -0.1225156446376204440720569753e+01,
    -0.4957589496572501915214079952e+00,
    0.1664377182454986536961530415e+01,
    -0.3503288487499736816886487290e+00,
    0.3341791187130174790297318841e+00,
    0.8192320648511571246570742613e-01,
    -0.2235530786388629525884427845e-01,
    0.1e+00,
    0.2e+00,
    0.777777777777777777777777777778e+00,
    -0.84289382761090128651353491142e+01,
    0.56671495351937776962531783590e+00,
    -0.30689499459498916912797304727e+01,
    0.23846676565120698287728149680e+01,
    0.21170345824450282767155149946e+01,
    -0.87139158377797299206789907490e+00,
    0.22404374302607882758541771650e+01,
    0.63157877876946881815570249290e+00,
    -0.88990336451333310820698117400e-01,
    0.18148505520854727256656404962e+02,
    -0.91946323924783554000451984436e+01,
    -0.44360363875948939664310572000e+01,
    0.10427508642579134603413151009e+02,
    0.24228349177525818288430175319e+03,
    0.16520045171727028198505394887e+03,
    -0.37454675472269020279518312152e+03,
    -0.22113666853125306036270938578e+02,
    0.77334326684722638389603898808e+01,
    -0.30674084731089398182061213626e+02,
    -0.93321305264302278729567221706e+01,
    0.15697238121770843886131091075e+02,
    -0.31139403219565177677282850411e+02,
    -0.93529243588444783865713862664e+01,
    0.35816841486394083752465898540e+02,
    0.19985053242002433820987653617e+02,
    -0.38703730874935176555105901742e+03,
    -0.18917813819516756882830838328e+03,
    0.52780815920542364900561016686e+03,
    -0.11573902539959630126141871134e+02,
    0.68812326946963000169666922661e+01,
    -0.10006050966910838403183860980e+01,
    0.77771377980534432092869265740e+00,
    -0.27782057523535084065932004339e+01,
    -0.60196695231264120758267380846e+02,
    0.84320405506677161018159903784e+02,
    0.11992291136182789328035130030e+02,
    -0.25693933462703749003312586129e+02,
    -0.15418974869023643374053993627e+03,
    -0.23152937917604549567536039109e+03,
    0.35763911791061412378285349910e+03,
    0.93405324183624310003907691704e+02,
    -0.37458323136451633156875139351e+02,
    0.10409964950896230045147246184e+03,
    0.29840293426660503123344363579e+02,
    -0.43533456590011143754432175058e+02,
    0.96324553959188282948394950600e+02,
    -0.39177261675615439165231486172e+02,
    -0.14972683625798562581422125276e+03,
    5.61675022830479523392909219681e-2,
    2.53500210216624811088794765333e-1,
    -2.46239037470802489917441475441e-1,
    -1.24191423263816360469010140626e-1,
    1.5329179827876569731206322685e-1,
    8.20105229563468988491666602057e-3,
    7.56789766054569976138603589584e-3,
    -8.298e-3,
    3.18346481635021405060768473261e-2,
    2.83009096723667755288322961402e-2,
    5.35419883074385676223797384372e-2,
    -5.49237485713909884646569340306e-2,
    -1.08347328697249322858509316994e-4,
    3.82571090835658412954920192323e-4,
    -3.40465008687404560802977114492e-4,
    1.41312443674632500278074618366e-1,
    -4.28896301583791923408573538692e-1,
    -4.69762141536116384314449447206e0,
    7.68342119606259904184240953878e0,
    4.06898981839711007970213554331e0,
    3.56727187455281109270669543021e-1,
    -1.39902416515901462129418009734e-3,
    2.9475147891527723389556272149e0,
    -9.15095847217987001081870187138e0)

abstract ODEProblem

type ODEProblemFunction
    f::Function
end

F!(p::ODEProblemFunction, y, x, t) = copy!(y, p.f(x, t))

dop853(p, y0, tspan; args...) = odedop(dop853coeff, p, y0, tspan; args...)
dopri5(p, y0, tspan; args...) = odedop(dopri5coeff, p, y0, tspan; args...)

dop853(f::Function, y0::Vector, tspan; args...) = dop853(ODEProblemFunction(f), y0, tspan; args...)
dopri5(f::Function, y0::Vector, tspan; args...) = dopri5(ODEProblemFunction(f), y0, tspan; args...)

function dop853(f::Function, y0::Number, tspan; args...)
    tout, yout = dop853(ODEProblemFunction(f), [y0], tspan; args...)
    return tout, vcat(yout...)
end

function dopri5(f::Function, y0::Number, tspan; args...)
    tout, yout = dopri5(ODEProblemFunction(f), [y0], tspan; args...)
    return tout, vcat(yout...)
end

function odedop(coeff, p, y0, tspan;
    reltol::Vector{Float64}=[1e-6], abstol::Vector{Float64}=[sqrt(eps())],
    safe::Float64=0.9, fac1::Float64=getfac1(coeff),
    fac2::Float64=getfac2(coeff), beta::Float64=getbeta(coeff),
    maxstep::Float64=tspan[end]-tspan[1], initstep=0.0,
    maxsteps::Int64=100000, printmessages::Bool=false, nstiff::Int64=1000,
    iout::Int64=0, solout::Function=s(x...)=return, dense::Vector{Int64}=[1:length(y0)],
    points::Symbol=:all)

    irtrn = 0

    x = tspan[1]
    xend = tspan[end]
    tout = Array(typeof(tspan[1]),0)
    yout = Array(typeof(y0*one(x)),0)

    if points == :all
        push!(tout, x)
        push!(yout, y0)
    end

    n = length(y0)
    nrd = length(dense)
    nfcn = 0
    nstep = 0
    naccpt = 0
    nrejct = 0
    k1 = 0.*y0
    k2 = 0.*y0
    k3 = 0.*y0
    k4 = 0.*y0
    k5 = 0.*y0
    k6 = 0.*y0
    k7 = 0.*y0
    k8 = 0.*y0
    k9 = 0.*y0
    k10 = 0.*y0
    y = 0.*y0

    copy!(y, y0)
    facold = 1e-4
    expo1 = getexpo1(coeff, beta)
    facc1 = 1.0/fac1
    facc2 = 1.0/fac2
    posneg = sign(xend-x)
    abstol = length(abstol) == 1 ? fill(abstol[1], n) : abstol
    reltol = length(reltol) == 1 ? fill(reltol[1], n) : reltol
    last = false
    hlamb = 0.0
    iasti = 0
    F!(p, k1, x, y)
    nfcn += 1
    hmax = abs(maxstep)
    nmax = maxsteps
    h = initstep
    iord = order(coeff)
    if h == 0.0
        h = hinit(n, p, x, y, xend, posneg, k1, k2, k3, iord, hmax, abstol, reltol)
        if printmessages
            println("hinit = $h")
        end
        nfcn += 1
    end
    reject = false
    xold = x
    if iout != 0
        #hout = 1.0
        irtrn = solout(naccpt+1, xold, x, y, con, icomp)
        if irtrn < 0
            return y
        end
    end

    while true
        if nstep > nmax
            error("Exit at x=$x. More than nmax=$nmax steps needed.")
        end
        if 0.1*abs(h) <= abs(x)*eps()
            error("Exit at x=$x. Step size too small, h=$h.")
        end
        if (x+1.01*h-xend)*posneg > 0.0
            h = xend-x
            last = true
        end
        nstep += 1
        if irtrn >= 2
            F!(p, k1, x, y)
        end

        y1, err, ncore = dopcore(coeff, n, p, x, y, h, k1, k2, k3, k4, k5, k6, k7,
        k8, k9, k10, abstol, reltol)
        xph = x+h
        nfcn += ncore

        fac11 = err^expo1
        fac = fac11/facold^beta
        fac = max(facc2, min(facc1, fac/safe))
        hnew = h/fac
        if err <= 1.0
            facold = max(err, 1e-4)
            naccpt += 1
            # Stiffness detection
            if mod(naccpt, nstiff) == 0 || iasti > 0
                nonsti = 0
                stnum = 0.0
                stden = 0.0
                for i = 1:n
                    if typeof(coeff) == DOP853
                        stnum += (k4[i] - k3[i])^2
                        stden += (k5[i] - y1[i])^2
                    elseif typeof(coeff) == DOPRI5
                        stnum += (k2[i] - k6[i])^2
                        stden += (y1[i] - k7[i])^2
                    end
                end
                if stden > 0.0
                    hlamb = abs(h)*sqrt(stnum/stden)
                end
                if hlamb > 6.1
                    nonsti = 0
                    iasti += 1
                    if iasti == 15
                        if printmessages
                            println("The problem seems to become stiff at x=$x")
                        end
                    end
                else
                    nonsti += 1
                    if nonsti == 6
                        iasti = 0
                    end
                end
            end
            # Final preparation for dense output
            event = iout == 3 && xout <= xph
            if iout == 2 || event
                for (i,j) in enumerate(dense)
                end
            end
            copy!(k1, k4)
            copy!(y, k5)
            xold = x
            x = xph
            if points == :all
                push!(tout, x)
                push!(yout, copy(k5))
            end
            # if
            # solout
            # end
            # Normal exit
            if last
                h = hnew
                idid = 1
                stats = ["function_calls"=>nfcn,
                "steps"=>nstep,
                "accepted"=>naccpt,
                "rejected"=>nrejct]
                if points == :last
                    return x, k5, stats
                end
                return tout, yout, stats
            end
            if abs(hnew) > hmax
                hnew = posneg*hmax
            end
            if reject
                hnew = posneg*min(abs(hnew),abs(h))
            end
            reject = false
        else
            hnew = h/min(facc1,fac11/safe)
            reject = true
            if naccpt >= 1
                nrejct += 1
            end
            last = false
        end
        h = hnew
    end
end

function denseout(ind, t, told, h, coeff, dense)
    if ~in(ind, dense)
        error("No dense output available for component: $ind.")
    else
        s = (t - told)/h
        s1 = 1.0 - s
    end
end

function hinit(n::Int64, p, x::Float64, y::Vector, xend::Float64, posneg::Float64, f0::AbstractVector, f1::AbstractVector, y0::AbstractVector, iord::Int64, hmax::Float64, abstol::Vector{Float64}, reltol::Vector{Float64})
    dnf = 0.0
    dny = 0.0
    for i = 1:n
        sk = abstol[i] + reltol[i]*abs(y[i])
        dnf += (f0[i]/sk)^2
        dny += (y[i]/sk)^2
    end
    if maximum(abs(dnf)) <= 1e-10 || maximum(abs(dny)) <= 1e-10
        h = 1e-6
    else
        h = sqrt(maximum(abs(dny/dnf)))*0.01
    end
    h = min(h, hmax)
    h = h*posneg
    y1 = y + h*y0
    F!(p, f1, x+h, y1)
    der2 = 0.0
    for i = 1:n
        sk = abstol[i] + reltol[i]*abs(y[i])
        der2 += ((f1[i]-f0[i])/sk)^2
    end
    der2 = sqrt(maximum(abs(der2)))/h
    der12 = max(maximum(abs(der2)), sqrt(maximum(abs(dnf))))
    if der12 <= 1e-15
        h1 = max(1e-6, abs(h)*1e-3)
    else
        h1 = (0.01/der12)^(1.0/iord)
    end
    h = min(100*abs(h), h1, hmax)
    return h*posneg
end

order(c::DOP853) = 8
order(c::DOPRI5) = 5
getfac1(c::DOP853) = 0.333
getfac1(c::DOPRI5) = 0.2
getfac2(c::DOP853) = 6.0
getfac2(c::DOPRI5) = 10.0
getbeta(c::DOP853) = 0.0
getbeta(c::DOPRI5) = 0.04
getexpo1(c::DOP853, beta) = 1.0/8.0 - beta*0.2
getexpo1(c::DOPRI5, beta) = 0.2 - beta*0.75

function dopcore(c::DOP853, n::Int64, p, x::Float64, y::Vector, h::Float64, k1::Vector, k2::Vector, k3::Vector, k4::Vector, k5::Vector, k6::Vector, k7::Vector, k8::Vector, k9::Vector, k10::Vector, abstol::Vector, reltol::Vector)
    y1 = 0.0 * y

    for i = 1:n
        y1[i] = y[i] + h*c.a21*k1[i]
    end
    F!(p, k2, x+c.c2*h, y1)
    for i = 1:n
        y1[i] = y[i] + h*(c.a31*k1[i] + c.a32*k2[i])
    end
    F!(p, k3, x+c.c3*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a41*k1[i]+c.a43*k3[i])
    end
    F!(p, k4, x+c.c4*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a51*k1[i]+c.a53*k3[i]+c.a54*k4[i])
    end
    F!(p, k5, x+c.c5*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a61*k1[i]+c.a64*k4[i]+c.a65*k5[i])
    end
    F!(p, k6, x+c.c6*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a71*k1[i]+c.a74*k4[i]+c.a75*k5[i]+c.a76*k6[i])
    end
    F!(p, k7, x+c.c7*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a81*k1[i]+c.a84*k4[i]+c.a85*k5[i]+c.a86*k6[i]+c.a87*k7[i])
    end
    F!(p, k8, x+c.c8*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a91*k1[i]+c.a94*k4[i]+c.a95*k5[i]+c.a96*k6[i]+c.a97*k7[i]+c.a98*k8[i])
    end
    F!(p, k9, x+c.c9*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a101*k1[i]+c.a104*k4[i]+c.a105*k5[i]+c.a106*k6[i]
           +c.a107*k7[i]+c.a108*k8[i]+c.a109*k9[i])
    end
    F!(p, k10, x+c.c10*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a111*k1[i]+c.a114*k4[i]+c.a115*k5[i]+c.a116*k6[i]
           +c.a117*k7[i]+c.a118*k8[i]+c.a119*k9[i]+c.a1110*k10[i])
    end
    F!(p, k2, x+c.c11*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a121*k1[i]+c.a124*k4[i]+c.a125*k5[i]+c.a126*k6[i]
           +c.a127*k7[i]+c.a128*k8[i]+c.a129*k9[i]+c.a1210*k10[i]+c.a1211*k2[i])
    end
    F!(p, k3, x+h, y1)
    for i = 1:n
        k4[i] = c.b1*k1[i]+c.b6*k6[i]+c.b7*k7[i]+c.b8*k8[i]+c.b9*k9[i]+c.b10*k10[i]+c.b11*k2[i]+c.b12*k3[i]
        k5[i] = y[i]+h*k4[i]
    end

    # Error estimation
    err = 0.0
    err2 = 0.0
    for i = 1:n
        sk = abstol[i] + reltol[i]*max(maximum(abs(y[i])), maximum(abs(k5[i])))
        erri = k4[i] - c.bhh1*k1[i] - c.bhh2*k9[i] - c.bhh3*k3[i]
        err2 += (erri/sk)*(erri/sk)
        erri = c.er1*k1[i] + c.er6*k6[i] + c.er7*k7[i] + c.er8*k8[i] + c.er9*k9[i] + c.er10*k10[i] + c.er11*k2[i] + c.er12*k3[i]
        err += (erri/sk)*(erri/sk)
    end
    deno = maximum(abs(err)) + 0.01*maximum(abs(err2))
    if deno <= 0.0
        deno = 1.0
    end
    err = abs(h)*err*sqrt(1.0/(n*deno))

    F!(p, k4, x+h, k5)

    # Number of function evaluations
    neval = 12
    return y1, maximum(abs(err)), neval
end

function dopcore(c::DOPRI5, n::Int64, p, x::Float64, y::Vector, h::Float64, k1::Vector, k2::Vector, k3::Vector, k4::Vector, k5::Vector, k6::Vector, k7::Vector, k8::Vector, k9::Vector, k10::Vector, abstol::Vector, reltol::Vector)
    y1 = 0.0 * y

    for i = 1:n
        y1[i] = y[i] + h*c.a21*k1[i]
    end
    F!(p, k2, x+c.c2*h, y1)
    for i = 1:n
        y1[i] = y[i] + h*(c.a31*k1[i] + c.a32*k2[i])
    end
    F!(p, k3, x+c.c3*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a41*k1[i]+c.a42*k2[i]+c.a43*k3[i])
    end
    F!(p, k8, x+c.c4*h, y1)
    for i = 1:n
        y1[i] = y[i]+h*(c.a51*k1[i]+c.a52*k2[i]+c.a53*k3[i]+c.a54*k8[i])
    end
    F!(p, k9, x+c.c5*h, y1)
    for i = 1:n
        k7[i] = y[i]+h*(c.a61*k1[i]+c.a62*k2[i]+c.a63*k3[i]+c.a64*k8[i]+c.a65*k9[i])
    end
    F!(p, k6, x+h, k7)
    for i = 1:n
        y1[i] = y[i]+h*(c.a71*k1[i]+c.a73*k3[i]+c.a74*k8[i]+c.a75*k9[i]+c.a76*k6[i])
    end
    F!(p, k4, x+h, y1)
    for i = 1:n
        k8[i] = (c.e1*k1[i]+c.e3*k3[i]+c.e4*k8[i]+c.e5*k9[i]+c.e6*k6[i]+c.e7*k4[i])*h
        k5[i] = y1[i]
    end

    # Error estimation
    err = 0.0
    for i = 1:n
        sk = abstol[i] + reltol[i]*max(maximum(abs(y[i])), maximum(abs(y1[i])))
        err += (k8[i]/sk)*(k8[i]/sk)
    end
    err = sqrt(err/n)

    # Number of function evaluations
    neval = 6
    return y1, maximum(abs(err)), neval
end
