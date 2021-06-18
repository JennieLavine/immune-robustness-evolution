using JSON

X = zeros(2, 3)
X[1,2] = true

println(JSON.json(X))
