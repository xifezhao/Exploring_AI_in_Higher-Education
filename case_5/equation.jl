# Define constants
ħ = 1.05457e-34  # Reduced Planck constant (J⋅s)
m = 9.109e-31  # Mass of electron (kg)
a = 1e-9        # Width of the well (m)

# Define potential energy function (infinite outside the well)
function V(x)
  if x < 0 || x > a
    return Inf
  else
    return 0.0
  end
end

# Define the Hamiltonian operator
function H(ψ)
  return -ħ^2 / (2*m) * diff(ψ, x)^2 + V(x) * ψ(x)
end

# Define boundary conditions (ψ(0) = 0, ψ'(a/4) = 0)
bc1(ψ) = ψ(0)
bc2(ψ) = diff(ψ, x)(a/4)

# Discretize the spatial domain
x = range(0.0, stop=a, length=100)  # 100 points between 0 and a

# Solve the eigenvalue problem with boundary conditions
eigenvals, eigenvecs = solveev(H, x, [bc1, bc2])

# Print the first few eigenvalues (energy levels)
println("First few eigenvalues (in Joules):")
println(eigenvals[1:5])
