#
# Aplica condições de contorno essenciais 
# homogêneas
# ! após o nome para identificar que está sendo alterado
# os vetores/matrizes fornecidos
function Aplica_CC_homo!(apoios,K,F) 

    # Aplica as condições de contorno essenciais homogêneas
    for i=1:size(apoios,1) 

        # No e gl local do apoio
        no  = Int(apoios[i,1])
        gll = Int(apoios[i,2])

        # Testa se o usuário tem ideia do que
        # ele está fazendo
        # Deslocamentos conhecidos precisam ser zero
        if apoios[i,3]!=0 
            error("Método de resolução incompátivel. Tente Multiplicadores de Lagrande ou Rearranjo")
        end

        # Gl global
        gl = 2*(no-1) + gll

        # Zera linha e coluna da rigidez
        K[gl,:] .= 0.0
        K[:,gl] .= 0.0

        # Zera o vetor de forças nesse gl
        F[gl] = 0.0

        # Coloca 1.0 na diagonal
        K[gl,gl] = 1.0

    end #gls

end


#
# Aplica cc essenciais por Multiplicadores de Lagrange
#
# Monta o sistema aumentado
# [K S' {U}= {F}
#  S 0] {L}  {Ub}
#
#
function Aplica_CC_lagrange(nnos,apoios,K,F)

    # Número de gls do problema original
    n = 2*nnos

    # Número de cc essenciais
    m = size(apoios,1)

    # Define o sistema aumentado de equações
    KA = zeros(n+m,n+m)
    FA = zeros(n+m)

    # Define a matriz S e o vetor Ub
    S  = zeros(m,n)
    Ub = zeros(m)

    # Aplica as condições de contorno essenciais
    for i=1:m 

        # No e gl local do apoio
        no  = Int(apoios[i,1])
        gll = Int(apoios[i,2])
        valor = apoios[i,3]

        # Gl global
        gl = 2*(no-1) + gll

        # Posiciona na linha da matriz S
        S[i,gl] = 1.0

        # Posiciona o valor em Ub
        Ub[i] = valor

    end

    # Posiciona os blocos na matriz e no vetor aumentados
    KA[1:n,1:n]     .= K
    KA[1:n,n+1:end] .= S'
    KA[n+1:end,1:n] .= S

    FA[1:n]     .= F
    FA[n+1:end] .= Ub

    return KA, FA
end

#
# Aplica condições de contorno essenciais (deslocamentos
# prescritos) e reslve por rearranjo 
#
function Rearranjo!(nnos,apoios,K,F)
# Numero de graus de liberdade do problema

    # graus de liberdado
    n = 2*nnos

    # Modifica o vetor de carregamentos
    for i=1:n

        for j=1:size(apoios,1)

        # vetor deslocamentos
        no  = Int(apoios[j,1])
        gll = Int(apoios[j,2])
        valor = Int(apoios[j,3])

        # Gl global
        gl = 2*(no-1) + gll
        
        # Subtrai a contribuição de cada apoio
        F[i] = F[i] - K[i,gl] * valor

        end #j
    end #i

    # Modifica matriz de rigidez zerando linha e co-
    # luna

    for i=1:size(apoios,1) 

        # No e gl local do apoio
        no  = Int(apoios[i,1])
        gll = Int(apoios[i,2])

        # Gl global
        gl = 2*(no-1) + gll

        # Zera linha e coluna da rigidez
        K[gl,:] .= 0.0
        K[:,gl] .= 0.0

        # Zera o vetor de forças nesse gl
        F[gl] = 0.0

        # Coloca 1.0 na diagonal
        K[gl,gl] = 1.0

    end #i

    return K, F
end

