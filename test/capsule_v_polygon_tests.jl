
function create_n_sided(N,d)
    ns = [ [cos(θ);sin(θ)] for θ = 0:(2*π/N):(2*π*(N-1)/N)]
    #
    A = vcat(transpose.((ns))...)
    b = d*ones(N)
    return SMatrix{N,2}(A), SVector{N}(b)
end

for i = 1:1000

    cap = DP.create_capsule(:MRP)
    cap.R = 0.7
    cap.L = 1.3

    A2, b2 = create_n_sided(6,1.2)
    R2 = abs(2*randn())
    polyg = DP.create_polygon(A2,b2,R2,:MRP)

    r1 = 2*(@SVector randn(3))
    r2 = 2*(@SVector randn(3))
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)
    cap.r = r1
    polyg.r = r2
    cap.p = p1
    polyg.p = p2

    if DP.proximity(cap,polyg) > 0
        p1,p2 = DP.closest_points(cap,polyg)
        @test norm(p1 - p2) > 0
    end

end
