# some standard test cases
const test_cases =
    Dict(:constant_in_time=>
         Dict(:ivp      => ExplicitODE(0.0,[0.0],
                                       (t,y,dy)->dy[:]=6.0,
                                       J! = (t,y,dy)->dy[1]=0.0),
              :sol      => t->[6t],
              :name     => "y'=6 (vector)",
              :options  => Dict(:initstep => 0.1,
                                :tstop => 1.0,
                                :tout  => [0.0,0.1,1.0])
              ),

         :variable_in_time=>
         Dict(:ivp    => ExplicitODE(0.0,[0.0],
                                     (t,y,dy)->dy[1]=2t,
                                     J! = (t,y,dy)->dy[1]=0.0),
              :sol    => t->[t^2],
              :name   => "y'=2t",
              :options=> Dict(:tout => [0:0.001:1;],
                              :tstop => 1.0,
                              :initstep => 0.001)
              ),

         :linear=>
         Dict(:ivp    => ExplicitODE(0.0,[1.0],
                                     (t,y,dy)->dy[1]=y[1],
                                     J! = (t,y,dy)->dy[1]=1.0),
              :sol    => t->[exp(t)],
              :name   => "y'=y",
              :options=> Dict(:tout => [0:0.001:1;],
                              :tstop => 1.0,
                              :initstep => 0.001)
              ),

         :backward_in_time=>
         Dict(:ivp    => ExplicitODE(1.0,[1.0],
                                     (t,y,dy)->dy[1]=y[1],
                                     J! = (t,y,dy)->dy[1]=1.0),
              :sol    => t->[exp(t-1)],
              :name   => "y'=y backwards",
              :options=> Dict(:tout => [1:-0.001:0;],
                              :tstop => 0.0,
                              :initstep => 0.001)
              ),

         :harmonic=>
         Dict(:ivp    => ExplicitODE(0.0,[1.0,2.0],
                                     (t,y,dy)->(dy[1]=-y[2];dy[2]=y[1]),
                                     J! = (t,y,dy)->copy!(dy,Float64[[0,1] [-1,0]])),
              :sol    => t->[cos(t)-2*sin(t), 2*cos(t)+sin(t)],
              :name   => "harmonic (with Jacobian)",
              :options=> Dict(:tout => [0:.1:1;],
                              :tstop => 1.0,
                              :initstep => 0.001)
              ),

         :harmonic_no_jac=>
         Dict(:ivp    => ExplicitODE(0.0,[1.0,2.0],
                                     (t,y,dy)->(dy[1]=-y[2];dy[2]=y[1])),
              :sol    => t->[cos(t)-2*sin(t), 2*cos(t)+sin(t)],
              :name   => "harmonic (no Jacobian)",
              :options=> Dict(:tout => [0:.1:1;],
                              :tstop => 1.0,
                              :initstep => 0.001)
              ),

         :harmonic_minimal_types=>
         Dict(:ivp      => ExplicitODE(MyFloat(0.0),
                                       Position(MyFloat(0.0),MyFloat(1.0)),
                                       (t,y,dy)->(dy.x=y.y;dy.y=-y.x),
                                       J! = (t,y,J)->(J[:]=0;J[2,1]=1;J[1,2]=-1)),
              :sol      => t->Position(sin(t),cos(t)),
              :name     => "harmonic (minimal types)",
              :options  => Dict(:initstep => MyFloat(0.001),
                                :tstop => MyFloat(1.0),
                                :tout  => [MyFloat(0.0),MyFloat(0.1),MyFloat(1.0)])
              ),

         # TODO: ForwardDiff Jacobian doesn't seem to work with custom AbstractVector type

         # :harmonic_minimal_types_no_jac=>
         # Dict(:ivp      => ExplicitODE(MyFloat(0.0),
         #                               Position(MyFloat(0.0),MyFloat(1.0)),
         #                               (t,y,dy)->(dy.x=y.y;dy.y=-y.x)),
         #      :sol      => t->Position(sin(t),cos(t)),
         #      :name     => "harmonic (minimal types, no Jacobian)",
         #      :options  => Dict(:initstep => MyFloat(0.1),
         #                        :tstop => MyFloat(1.0),
         #                        :tout  => [MyFloat(0.0),MyFloat(0.1),MyFloat(1.0)])
         #      )

         )
