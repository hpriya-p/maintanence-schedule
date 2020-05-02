FROM julia:latest
COPY src .
COPY instances .

# Add Julia libraries
ADD REQUIRE /.julia/v1.1/REQUIRE
RUN julia -e "using Pkg; Pkg.resolve()"
RUN julia solve_ip.jl "example1.json"