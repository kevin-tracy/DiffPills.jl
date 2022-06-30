# DiffPills.jl

<p align="center">
  <img width="600" src="https://github.com/kevin-tracy/DiffPills.jl/blob/master/extras/images/diff_pills_logo.png">
</p>

DiffPills.jl is a tool for computing proximity between two pill-shaped objects, and providing analytical smooth gradients of this operation.

## basic usage

```julia
import DifferentialProximity as dp

# create the pill objects (with :MRP or :quat for attitude)
pill_1 = create_pill(:MRP)
pill_2 = create_pill(:MRP)

# put position, attitude (MRP), radius, and length in
pill_1.r = @SVector randn(3);   pill_2.r = @SVector randn(3)
pill_1.p = @SVector randn(3);   pill_2.p = @SVector randn(3)
pill_1.R = 0.75;                pill_1.R = 1.5
pill_1.L = 3.0;                 pill_1.L = 2.0

# calculate proximity
ℓ = dp.pill_proximity(pill_1, pill_2)

# calculate proximity gradients
∂ℓ_∂r1, ∂ℓ_∂p1, ∂ℓ_∂r2, ∂ℓ_∂p2 = dp.pill_proximity_jacobians(pill_1, pill_2)
```
