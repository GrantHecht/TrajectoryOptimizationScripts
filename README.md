# TrajectoryOptimizationScripts

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> TrajectoryOptimizationScripts

It is authored by Grant Hecht.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Download the dependancies which are not included in the Julia package repository,
   i.e., GrantHecht/AstroUtils and GrantHecht/AstroEOMs. 
   These can be placed anywhere, but I would recommend the file path:
   /home/user-name/.julia/dev/AstroUtils and
   /home/user-name/.julia/dev/AstroEOMs
1. Open a Julia console and do:
   ```
   julia> ] # enter package mode
   (@v1.8) pkg> add DrWatson # install globally, for using `quickactivate`
   (@v1.8) pkg> activate path/to/this/project
   (@v1.8) pkg> dev /home/user-name/.julia/dev/AstroUtils /home/user-name/.julia/dev/AstroEOMs
   (@v1.8) pkg> instantiate # afterwards, hit backspace to exit package mode
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "TrajectoryOptimizationScripts"
```
which auto-activate the project and enable local path handling from DrWatson.
