FROM julia:latest
COPY src .
COPY instances .
COPY REQUIRE .
# Add Julia libraries
RUN for i in $(cat REQUIRE); do julia -e "using Pkg; Pkg.add($i); Pkg.build($i)"; done

RUN julia -e "using Pkg; Pkg.resolve()"
RUN julia solve_ip.jl "example1.json"