


function calc_W_from_B(Tk, B, P)

    pws = Pws(B)
    efac = efactor(B,P)
    xsv = efac * pws / P
    w2 = Mv/Ma * xsv/(1-xsv)

    w = (enthalpyair(B,P) - enthalpyair(Tk,P) - w2*(enthalpywi(B) - enthalpyvapor(B))) /
        (enthalpyvapor(Tk) - enthalpywi(B))
    
    f0(w0) = aux_WB(w0, Tk, B, P, pws, efac)
    return Roots.find_zero(f0,w)
end

function aux_WB(w, Tk, B, P, pws, efac)
    xv1 = w / (Mv/Ma+w)
    xv2 = efac * pws / P
    w2 = Mv / Ma * xv2 / (1.0-xv2)
    #(1.0+w)*enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - (1.0+w2)*enthalpymoist(B,P,xv2)
    enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - enthalpymoist(B,P,xv2)
end

function aux_WB(w, Tk, B, P)
    pws = Pws(B)
    efac = efactor(B,P) 
    xv1 = w / (Mv/Ma+w)
    xv2 = efac * pws / P
    w2 = Mv / Ma * xv2 / (1.0-xv2)
    #(1.0+w)*enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - (1.0+w2)*enthalpymoist(B,P,xv2)
    enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - enthalpymoist(B,P,xv2)
end

function calcwetbulb(Tk, P, xv)
    w = humrat(xv)
    B = Tk - 1.0 # Initial guess
    h = 1e-7
    dB = 0.0
    f0(_B) = aux_WB(w, Tk,_B, P)
    return Roots.find_zero(f0,B)     
end

humrat(xv) = xv / (1-xv) * (Mv/Ma)