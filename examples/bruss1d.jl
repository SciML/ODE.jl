using ODE,DASSL

const T = Float64 # the datatype used, probably Float64
Tarr = SparseMatrixCSC
const N = 500 # active gridpoints
const dof = 2N # degrees of freedom (boundary points are not free)
dx = 1/(N+1)
x = dx:dx:1-dx
dae = 0  # index of DAE, ==0 for ODE

# parameters
alpha = 1/50
const gamma = alpha/dx^2

# BC
const one_ = one(T)
const three = 3*one_
ubc1() = one_
ubc2() = one_
vbc1() = three
vbc2() = three
#IC
u0(x) = 1 + 0.5*sin(2π*x)
v0(x) = 3 * ones(length(x))

inds(i) = (2i-1, 2i)
function fn!(res, y, dydt)
    # the ode function res = f(y,dydt)
    #
    # Note the u and v are staggered:
    # u = y[1:2:end], v=y[2:2:end]
    
    # setup for first step
    i = 1
    iu, iv = inds(i)
    ui0 = ubc1(); vi0 = vbc1() # BC
    ui1 = y[iu];   vi1 = y[iv]
    ui2 = y[iu+2]; vi2 = y[iv+2]
    while i<=N
        res[iu] = 1 + ui1^2*vi1 - 4*ui1 + gamma*(ui0 - 2*ui1 + ui2) - dydt[iu]
        res[iv] = 3*ui1 - ui1^2*vi1     + gamma*(vi0 - 2*vi1 + vi2) - dydt[iv]
        # setup for next step:
        i += 1
        iu, iv = inds(i)
        ui0 = ui1; vi0 = vi1
        ui1 = ui2; vi1 = vi2
        if i<N
            ui2 = y[iu+2]
            vi2 = y[iv+2]
        else # on last step do BC
            ui2 = ubc2()
            vi2 = vbc2()
        end
    end
    return nothing
end
fn!(;T_::Type=T, dof_=dof) = zeros(T_,dof_)
function fn(t, y, ydot)
    res = fn!()
    fn!(res, y, ydot)
    return res
end
function fn_ex(t, y)
    res = fn!()
    fn!(res, y, 0*y)
    return res
end



# Jacobian is a banded matrix with upper and lower bandwidth 2,
# c.f. W&H p.148, but Julia does not support banded matrices yet.
# Thus use a sparse matrix.

function jac!(dfdy, y, ydot, a)
    # The Jacobian of fn: = df/dy + a* df/dydot
    ii = 1 # direct index into dfdy.nzval
    for i=1:N # iterate over double-columns
        iu, iv = inds(i)
        ui = y[iu]
        vi = y[iv]
        # iu: column belonging to d/du(i)
        if iu-2>=1
            # du'(i-1)/du(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        if iv-2>=1
            # dv'(i-1)/du(i)==0
            dfdy.nzval[ii] = 0
            ii +=1
        end
        # du'(i)/du(i)
        dfdy.nzval[ii]   = 2*ui*vi - 4 -2*gamma - a
        ii += 1
        # dv'(i)/du(i)
        dfdy.nzval[ii] = 3 - 2*ui*vi
        ii +=1
        if  iu+2<=2N
            # du'(i+1)/du(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        
        # iv: column belonging to d/dv(i)
        if iv-2>=1
            # dv'(i-1)/dv(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        # du'(i)/dv(i)
        dfdy.nzval[ii] = ui^2
        ii +=1
        # dv'(i)/dv(i)
        dfdy.nzval[ii]   = -ui^2 -2*gamma - a
        ii +=1
        if  iu+2<=2N
            # du'(i+1)/dv(i)==0
            dfdy.nzval[ii] = 0
            ii +=1
        end
        if iv+2<=2N
            # dv'(i+1)/dv(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
    end
    return nothing
end
function jac!(;T_::Type=T, dof_=dof)
    B = (ones(T_,dof_-2), ones(T_,dof_-1), ones(T_,dof_), ones(T_,dof_-1), ones(T_,dof_-2))        
    spdiagm(B, (-2,-1,0,1,2))
end
function jac(t, y, ydot, a)
    J = jac!()
    jac!(J, y, ydot, a)
    return J
end
function jac_ex(t, y)
    return jac(t, y, y, 0)
end

refsol = T[0.9949197002317599, 3.0213845767604077, 0.9594350193986054, 3.0585989778165419, 0.9243010095428502, 3.0952478919989637, 0.8897959106772672,
           3.1310118289054731, 0.8561653620284367, 3.1656101198770159, 0.8236197147449046, 3.1988043370624344, 0.7923328094811884, 3.2303999530641514,
           0.7624421042573115, 3.2602463873623941, 0.7340499750795348, 3.2882356529108807, 0.7072259700779899, 3.3142998590079271, 0.6820097782458483,
           3.3384078449410937, 0.6584146743834650, 3.3605612157873943, 0.6364312187752559, 3.3807900316323134, 0.6160310186921587, 3.3991483695914764,
           0.5971703941198909, 3.4157099395342736, 0.5797938277687891, 3.4305638938070224, 0.5638371159206763, 3.4438109320334580, 0.5492301695479158,
           3.4555597666485198, 0.5358994429426996, 3.4659239846027008, 0.5237699892215797, 3.4750193162238476, 0.5127671585747183, 3.4829613034792271,
           0.5028179665048467, 3.4898633463634923, 0.4938521662914935, 3.4958350971335204, 0.4858030633656755, 3.5009811668111510, 0.4786081100251151,
           3.5054001059792705, 0.4722093177200750, 3.5091836216744015, 0.4665535216425440, 3.5124159935026285, 0.4615925290790646, 3.5151736544621075,
           0.4572831793403656, 3.5175249049438184, 0.4535873393501199, 3.5195297317024448, 0.4504718553589467, 3.5212397070273984, 0.4479084778719241,
           3.5226979467564341, 0.4458737738041973, 3.5239391090719634, 0.4443490371324889, 3.5249894191569453, 0.4433202068820853, 3.5258667077466495,
           0.4427777991494095, 3.5265804544017270, 0.4427168579654424, 3.5271318289682063, 0.4431369281018266, 3.5275137272135266, 0.4440420513508381,
           3.5277107990730161, 0.4454407863109616, 3.5276994703501980, 0.4473462502188303, 3.5274479611304068, 0.4497761798232572, 3.5269163066324394,
           0.4527530066369863, 3.5260563887768472, 0.4563039400688689, 3.5248119894251024, 0.4604610498812091, 3.5231188790654930, 0.4652613370907894,
           3.5209049576992761, 0.4707467798082714, 3.5180904678044698, 0.4769643375804777, 3.5145883024867057, 0.4839658945842979, 3.5103044351908528,
           0.4918081185812277, 3.5051385005173827, 0.5005522089940899, 3.4989845585737802, 0.5102635039989190, 3.4917320776245013, 0.5210109134090777,
           3.4832671712209993, 0.5328661417420937, 3.4734741260299615, 0.5459026646938675, 3.4622372546582585, 0.5601944229089820, 3.4494431032230182,
           0.5758142001453760, 3.4349830354873361, 0.5928316594749734, 3.4187562033108012, 0.6113110218368440, 3.4006728962523969, 0.6313083867734524,
           3.3806582409098729, 0.6528687160104193, 3.3586561928427350, 0.6760225267555723, 3.3346337311179157, 0.7007823726569721, 3.3085851288057930,
           0.7271392249346637, 3.2805361342349380, 0.7550589020044152, 3.2505478606008622, 0.7844787296769868, 3.2187201496972175, 0.8153046416214843,
           3.1851941538893653, 0.8474089465959840, 3.1501538739882800, 0.8806289904192589, 3.1138264039027113, 0.9147669230929857, 3.0764806689389470,
           0.9495907429372025, 3.0384245041548366, 0.9848367306701233] # reference sol from Hairer has 1.36 scd
refsolinds = collect(1:7:dof) # refsol is only at components 1:7:2N


ic = zeros(T,2N)
ic[1:2:2N] = u0(x)
ic[2:2:2N] = v0(x)
tspan = T[0.0, 10.0] # integration interval
tspan0 = T[0.0, 0.001] # integration interval

#t,ydassl = dasslSolve(fn, ic, tspan, jacobian=jac)

## ode23s
# tout, yout = ode23s(fn_ex, ic, tspan0, jacobian=jac_ex)
# @time tout, yout = ode23s(fn_ex, ic, tspan, jacobian=jac_ex)
# @show norm(abs(yout[end][refsolinds]-refsol)./refsol, Inf)

# ## DASSL
# t,ydassl = dasslSolve(fn, ic, tspan0, jacobian=jac)
# @time t,ydassl = dasslSolve(fn, ic, tspan, jacobian=jac)
# @show norm(abs(ydassl[end][refsolinds]-refsol)./refsol, Inf)

## ROSW

# # No jac works but is slower and gives much higher precision than
# # other methods.
# tout, yout = ode_rosw(fn!, ic, tspan0)
# @time tout, yout = ode_rosw(fn!, ic, tspan)
# @show norm(abs(yout[end][refsolinds]-refsol)./refsol, Inf)

# Analytic Jacobian
tout, yout = ode_rosw(fn!, ic, tspan0, jacobian=(jac!, jac!()))
@time tout, yout = ode_rosw(fn!, ic, tspan, jacobian=(jac!, jac!()))
@show norm(abs(yout[end][refsolinds]-refsol)./refsol, Inf)

# Numerical colored Jacobian
tout, yout = ode_rosw(fn!, ic, tspan0, jacobian=jac!())
@time tout, yout = ode_rosw(fn!, ic, tspan, jacobian=jac!())
@show norm(abs(yout[end][refsolinds]-refsol)./refsol, Inf)

