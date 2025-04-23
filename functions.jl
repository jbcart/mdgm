import StatsBase.sample
import StatsBase.Weights

# function to create neighborhood structure for a square lattice
function generatelattice(n::Integer, dependence::Integer=1)
  if !(dependence ∈ [1, 2])
    throw(ArgumentError("dependence must be 1 or 2, corresponding to 1st or 2nd order dependence"))
  end
  if n < 1
    throw(ArgumentError("n must be >= 1"))
  end
  N = n^2 # total number of sites on the lattice
  lattice = [Vector{Integer}(undef, 1) for _ in 1:N]
  if dependence == 2
    temporary = zeros(Integer, 8)
  else # if dependence == 1
    temporary = zeros(Integer, 4)
  end
  for i in 1:N
    # check is i is (true=) not on ______ boundary
    l = rem(i, n) != 1             # left
    r = rem(i, n) != 0             # right
    t = div(i - 1, n) != (n - 1)   # top
    b = div(i - 1, n) != 0         # bottom
    i_t = i + n # index above i
    i_b = i - n # index below i
    if dependence == 2
      temporary .= [(i_b - 1) * b * l, i_b * b, (i_b + 1) * b * r,
										(i - 1) * l, (i + 1) * r,
										(i_t - 1) * t * l, i_t * t, (i_t + 1) * t * r]
    else # if dependence == 1
      temporary .= [i_b * b, (i - 1) * l, (i + 1) * r, i_t * t]
    end
    lattice[i] = temporary[temporary.!=0]
  end
  lattice
end

# Coupling from the past in Julia
# auxilary update_config_exact function
# α - scalar
function update_config_exact!(config, w, α::Number, ψ, ug, u, colors, coding)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      for k in 1:2
        w[k] = exp(α * colors[k])
      end
    elseif coding ∈ ["ising", "potts"]
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α * colors[k] + ψ * s)
      end
    elseif coding == "centered-auto"
      w[1] = 1 
      s = 0.0
      for j in eachindex(neighbors)
        s += neighbors[j] - exp(α)/(1+exp(α))
      end
       w[2] = exp(α + ψ * s)
    else 
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] * neighbors[j]
        end
        w[k] = exp(α * colors[k] + ψ * s)
      end
    end
    config[i] = u[i] < w[1] / sum(w) ? colors[1] : colors[2]
  end
end

# α - vector
function update_config_exact!(config, w, α::Vector, ψ, ug, u, colors, coding)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      for k in 1:2
        w[k] = exp(α[k] * colors[k])
      end
    elseif coding ∈ ["ising", "potts"]
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α[k] * colors[k] + ψ * s)
      end
    else 
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] * neighbors[j]
        end
        w[k] = exp(α[k] * colors[k] + ψ * s)
      end
    end
    config[i] = u[i] < w[1] / sum(w) ? colors[1] : colors[2]
  end
end

# α - Matrix
function update_config_exact!(config, w, α::Matrix, ψ, ug, u, colors, coding)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      for k in 1:2
        w[k] = exp(α[i,k] * colors[k])
      end
    elseif coding ∈ ["ising", "potts"]
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α[i,k] * colors[k] + ψ * s)
      end
    elseif coding == "centered-auto"
      w[1] = 1 
      s = 0.0
      for j in eachindex(neighbors)
        s += neighbors[j] - exp(α[ug[i][j],2])/(1+exp(α[ug[i][j],2]))
      end
      w[2] = exp(α[i,2] + ψ * s)
    else 
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] * neighbors[j]
        end
        w[k] = exp(α[i,k] * colors[k] + ψ * s)
      end
    end
    config[i] = u[i] < w[1] / sum(w) ? colors[1] : colors[2]
  end
end

# Gibbs update
function update_config_gibbs!(config, α::Number, ψ, ug, colors, coding)
  w = zeros(2)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      for k in 1:2
        w[k] = exp(α * colors[k])
      end
    elseif coding ∈ ["ising", "potts"]
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α * colors[k] + ψ * s)
      end
    elseif coding == "centered-auto"
      w[1] = 1 
      s = 0.0
      for j in eachindex(neighbors)
        s += neighbors[j] - exp(α)/(1+exp(α))
      end
      w[2] = exp(α + ψ * s)
    else 
      for k in 1:2
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] * neighbors[j]
        end
        w[k] = exp(α * colors[k] + ψ * s)
      end
    end
    config[i] = rand() < w[1] / sum(w) ? colors[1] : colors[2]
  end
end

# α - vector 
function update_config_gibbs!(config, α::Vector, ψ, ug, colors, coding)
  if coding ∈ ["autologistic", "ising-physics"]
    throw(ArgumentError("when α is a vector, coding must be ising or potts"))
  end
  q = length(colors)
  w = zeros(q)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      for k in 1:q
        w[k] = exp(α[k])
      end
    else
      for k in 1:q
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α[k] + ψ * s)
      end
    end
    config[i] = sample(colors, Weights(w))
  end
end

# α - Matrix
function update_config_gibbs!(config, α::Matrix, ψ, ug, colors, coding)
  q = length(colors)
  w = zeros(q)
  for i in eachindex(config)
    neighbors = @view config[ug[i]]
    if length(neighbors) == 0
      if coding ∈ ["autologistic", "ising-physics"]
        for k in 1:q
          w[k] = exp(α[i,k]*colors[k])
        end
      else
        for k in 1:q
          w[k] = exp(α[i,k])
        end
      end
    elseif coding ∈ ["autologistic", "ising-physics"]
      for k in 1:q
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] * neighbors[j]
        end
        w[k] = exp(α[i,k] * colors[k] + ψ * s)
      end
    elseif coding == "centered-auto"
      w[1] = 1 
      s = 0.0
      for j in eachindex(neighbors)
        s += neighbors[j] - exp(α[ug[i][j],2])/(1+exp(α[ug[i][j],2]))
      end
      w[2] = exp(α[i,2] + ψ * s)
    else
      for k in 1:q
        s = 0.0
        for j in eachindex(neighbors)
          s += colors[k] == neighbors[j]
        end
        w[k] = exp(α[i,k] + ψ * s)
      end
    end
    config[i] = sample(colors, Weights(w))
  end
end

# Swendsen-Wang Aux Functions
function update_clusters_sw!(clusters, visited, config, ug, ψ)
  c = 0
  for i in eachindex(config)
    if !visited[i]
      c += 1
      stack = [i]
      visited[i] = true

      while !isempty(stack)
        v = pop!(stack)
        clusters[v] = c
        for j in ug[v]
          if !visited[j] && config[j] == config[v] && rand() < (1-exp(-ψ))
            push!(stack,j)
            visited[j] = true
          end
        end
      end
    end
  end
  c
end

function update_config_sw!(config, clusters, visited, ug, α::Matrix, ψ, colors)
  q = length(colors)
  w = zeros(q)

  c = update_clusters_sw!(clusters, visited, config, ug, ψ)
  
  for i in 1:c
    for k in 1:q 
      w[k] = exp(sum(α[clusters .== i, k]))
    end
    if any(w .== Inf)
      for k in 1:q
        w[k] = w[k] == Inf ? 1 : 0
      end
    end
    config[clusters .== i] .= sample(colors, Weights(w))
  end
end

function update_config_sw!(config, clusters, visited, ug, α::Vector, ψ, colors)
  q = length(colors)
  w = zeros(q)
  for k in 1:q 
    w[k] = exp(α[k])
  end

  c = update_clusters_sw!(clusters, visited, config, ug, ψ)
  for i in 1:c
    config[clusters .== i] .= sample(colors, Weights(w))
  end
end

function update_config_sw!(config, clusters, visited, ug, α::Number, ψ, colors)
  w = exp.(α .* colors)

  c = update_clusters_sw!(clusters, visited, config, ug, ψ)
  
  for i in 1:c
    config[clusters .== i] .= sample(colors, Weights(w))
  end
end

function get_q(α::Matrix)
  size(α,2)
end

function get_q(α::Vector)
  length(α)
end

function get_q(α::Number)
  2
end

function sample_config(α::Number, colors, n)
  w = [1, exp(α)]
  sample(colors, Weights(w), n)
end

function sample_config(α::Vector, colors, n)
  w = exp.(α)
  sample(colors, Weights(w), n)
end

function sample_config(α::Matrix, colors, n)
  config = zeros(Int, n)
  for i in 1:n
    w = exp.(α[i,:])
    config[i] = sample(colors, Weights(w))
  end
  config
end

function rmrf(ug, α, ψ, M=50; coding="ising", method="exact", suppressmessage=false, maxsweeps=0)
  
  if method == "exact" && coding == "potts"
    throw(ArgumentError("Method='exact' only works for binary MRFs
      coding ∈ {'ising', 'autologistic', 'ising-physics'}"))
  end

  q = get_q(α)
  n = length(ug)

  if coding == "ising-physics"
    colors = [-1, 1]
  else 
    colors = 0:(q-1)
  end
	
  if method == "exact"
    w = zeros(2)
    initA = repeat([colors[1]], n)
    initB = repeat([colors[2]], n)
    configA = repeat([colors[1]], n)
    configB = repeat([colors[2]], n)
    sweeps = 1
    U = zeros(n,0)
    B = copy(M)
    while any(configA .!= configB)
      if maxsweeps != 0 && sweeps > maxsweeps
        println("max sweeps reached")
        break
      end
      if sweeps > 2
        M *= 2
      end
      if sweeps > 1
        B *= 2
      end
      if (!suppressmessage)
        println("Sweep $sweeps: $B iterations")
      end
      U = hcat(rand(n, M), U)
      configA .= initA
      configB .= initB
      for i ∈ 1:B
        u = @view U[:, i]
        update_config_exact!(configA, w, α, ψ, ug, u, colors, coding)
        update_config_exact!(configB, w, α, ψ, ug, u, colors, coding)
      end
      sweeps += 1
    end
    return configA
  else
    config = sample_config(α, colors, n)
    if method == "gibbs"
      for _ in 1:M
        update_config_gibbs!(config, α, ψ, ug, colors, coding)
      end
    elseif method == "sw"
      visited = repeat([false],n)
      clusters = zeros(n)
      for _ in 1:M
        update_config_sw!(config, clusters, visited, ug, α, ψ, colors)
        visited .= false
      end
    elseif method == "sw-gibbs"
      visited = repeat([false],n)
      clusters = zeros(n)
      for i in 1:M
        if i % 2 == 0
          update_config_sw!(config, clusters, visited, ug, α, ψ, colors)
          visited .= false
        else
          update_config_gibbs!(config, α, ψ, ug, colors, coding)
        end
      end
    end
    return config
  end
end

# function for pairwise sufficient statistic
function mrf_pairwise_ss(ug, config, coding="ising")
  if coding ∈ ["potts", "ising"]
    suffstat_i = map((x, i) -> sum(config[x] .== i), ug, config)
  else
    suffstat_i = map((x, i) -> sum(config[x] .* i), ug, config)
  end
  sum(suffstat_i) / 2
end

