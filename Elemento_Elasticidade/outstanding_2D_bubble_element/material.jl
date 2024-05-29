
# Monta a matriz constitutiva dado E,v e a hipótese
function Calcula_C(E,v,hipotese="EPT")

    # Testa se hipotese é válida
    hipotese in ["EPT","EPD"] || throw("Hipotese errada")

    G = E/(2*(1+v))
       
    if hipotese=="EPT"
        a = E/(1-v^2)
        C = [a v*a 0.0 ;
            v*a a  0.0 ;
              0  0  G]

    elseif hipotese=="EPD"
        d = (v+1)*(2*v-1)  
        a = E*(v-1)/d
        b = -(E*v)/d
        C = [a b 0 ;
             b a 0 ;
             0 0 G]
    else
 
        C = [a b b 0 0 0 ;
             b a b 0 0 0 ;
             b b a 0 0 0 ;
             0 0 0 G 0 0 ;
             0 0 0 0 G 0 ;
             0 0 0 0 0 G]

    end
    
    return C

end
