# roche_critical_points

Calculating and browsing over critical points of the Kopal potential Omega associated to a misaligned binary star.

Kopal potential is defined

    Omega(x,y,z, params) = 1/r1 + q(1/r2 - x/delta^2) + 1/2 (1 + q) F^2 [(x cos theta - z sin theta)^2 + y^2]

with

    r1 = sqrt(x^2 + y^2 + z^2)
    r2 = sqrt((x-delta)^2 + y^2 + z^2)

and the critical point r is determined by equation 
    
    nabla Omega(r) = 0;

Potential parameters are

    theta - angle of misalignment between orbit and star's angular momentum
    F - synchronicity parameter
    delta -  fractional instantaneous separation
    q - mass ratio

Content:

    /src 
        main.cpp        <- program for calculating critical points for a range of parameters 
        main_gpu.cu     <- a GPU version of the program
        main.h         
        Makefile

    /res
        plot.py         <- program for browsing over results
        res.pkl         <- compressed results
        bzip2_pickle.py


This software extens the content of the paper:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [Horvat et. al. "Physics of Eclipsing Binaries. III. Spin-Orbit Misalignment", ApJS 237:26 (14pp), 2018 August](https://doi.org/10.3847/1538-4365/aacd0f), preprint on [arxiv.org](https://arxiv.org/abs/1806.07680v2)

and is connected to [Project PHOEBE](http://phoebe-project.org/) and the GitHub repository https://github.com/phoebe-project/phoebe2.
