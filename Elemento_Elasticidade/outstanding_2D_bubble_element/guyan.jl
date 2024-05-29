#
# Aplica Condesação Estática
# (Redução de Guyan)
#
using LinearAlgebra

function Guyan(K::Matrix, glm::Vector{Int})
    # Numero total de graus de liberdade
    n = size(K,1)

    # Identifica os graus de liberdade
    gls = collect(1:n)

    # Identifica os graus de liberdade escravos
    gle = setdiff(gls, glm)

    # Dividindo a matriz K
    Kmm = K[glm, glm]
    Kme = K[glm, gle]
    Kem = K[gle, glm]
    Kee = K[gle, gle]

    # Calcula a matriz de rigidez reduzida
    # K = Kmm - Kme * inv(Kee) * Kem
    # Para fugir da inversão de matrix:
    X = Kee \ Kem
    K = Kmm - Kme*X

    return K

end