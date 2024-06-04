#
# conectividades é uma matriz com ne linhas (número de elementos)
# e 2 colunas (nó inicial e nó final)
#
function Global(ne,nnos,conectividades,coord, VE,Vnuxy, esp)

    # Precisamos definir a matriz global
    K = zeros(2*nnos,2*nnos)

    # Aloca vetores para as coordenadas de cada elemento
    X = zeros(8)
    Y = zeros(8)
    nos = zeros(8)
    gls = zeros(Int,16)

    # Loop nos elementos da malha
    for ele=1:ne

        # Recupera as informações do elemento
        Ee = VE[ele]
        ve = Vnuxy[ele]

        # Monta a matriz constitutiva do elemento
        C = Calcula_C(Ee,ve,"EPT")

        # Recupera as coordenadas do elemento
        for i=1:8
            # Descobre o nó
            no = conectividades[ele,i]
            
            # Guarda os nós do elemento
            nos[i] = no

            # Descobre as coordenadas
            X[i],Y[i] = coord[no,:]
        end

        # Monta a matriz de rigidez do elemento
        # no sistema local
        Kg = K_bi8n(C,X,Y,esp)

        # Agora precisamos posicionar Kg na matriz global do problema

        # Vetor com os gls GLOBAIS do elemento
        contador = 1
        # Varre os nós
        for i=1:8
            # Varre os gls do nó
            for j=1:2
                gls[contador] = Int(2*(nos[i]-1)+j)
                contador += 1
            end
        end

        # Soma Kg nas posições gls em K
        K[gls,gls] .= K[gls,gls] .+ Kg

    end #ele

    # Retorna a matriz de rigidez do problema
    return K

end


#
# Monta o vetor de força global
#
function Forca_global(nnos,forcas)

    # Montar o vetor de forças global
    F = zeros(2*nnos)

    # Para cada informação em forças 
    # posiciona a força no vetor de forças
    # globais do problema
    for i=1:size(forcas,1)
       
        # Recupera nó e gl local 
        no  = Int(forcas[i,1])
        gll = Int(forcas[i,2])

        # Recupera valor
        valor = forcas[i,3]

        # Adiciona ao valor da força
        F[2*(no-1)+gll] += valor

    end #i

    return F
end