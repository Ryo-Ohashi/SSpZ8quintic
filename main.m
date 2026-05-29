function truncated_gauss(K,a,b,c,d)
    poch_ab := K!1; poch_c1 := K!1; coef := [1];
    for n in [1..d] do
        poch_ab *:= (a+n-1)*(b+n-1); poch_c1 *:= (c+n-1)*n;
        coef := Append(coef,poch_ab/poch_c1);
    end for;
    return Polynomial(coef);
end function;

function SspZ8quintic(p)
    K := GF(p);
    R<x> := PolynomialRing(K);
    if p mod 8 eq 1 then
        g1 := truncated_gauss(K,1/2,1/4,3/4,(p-1)/4);
        g2 := truncated_gauss(K,1/4,1/8,7/8,(p-1)/8);
        g3 := truncated_gauss(K,3/4,1/8,3/8,(p-1)/8);
        g4 := truncated_gauss(K,3/4,3/8,5/8,(3*p-3)/8);
        g5 := truncated_gauss(K,3/4,9/8,11/8,(p-9)/8);
        g6 := truncated_gauss(K,3/4,11/8,13/8,(3*p-11)/8);
        F := [R!g1,g2,g3,g4,g5,g6];
    elif p mod 8 eq 3 then
        g1 := truncated_gauss(K,1/2,3/4,5/4,(p-3)/4);
        g2 := truncated_gauss(K,1/4,1/8,7/8,(3*p-1)/8);
        g3 := truncated_gauss(K,1/4,3/8,9/8,(p-3)/8);
        g4 := truncated_gauss(K,3/4,3/8,5/8,(p-3)/8);
        g5 := truncated_gauss(K,3/4,11/8,13/8,(p-11)/8);
        F := [R!g1,g2,g3,g4,g5];
    elif p mod 8 eq 5 then
        g1 := truncated_gauss(K,1/2,1/4,3/4,(p-1)/4);
        g2 := truncated_gauss(K,3/4,1/8,3/8,(5*p-1)/8);
        g3 := truncated_gauss(K,3/4,5/8,7/8,(p-5)/8);
        g4 := truncated_gauss(K,3/4,7/8,9/8,(3*p-7)/8);
        g5 := truncated_gauss(K,3/4,9/8,11/8,(5*p-9)/8);
        F := [R!g1,g2,g3,g4,g5];
    else
        g1 := truncated_gauss(K,1/2,3/4,5/4,(p-3)/4);
        g2 := truncated_gauss(K,3/4,7/8,9/8,(p-7)/8);
        F := [R!g1,g2];
    end if;
    G := GCD(F);
    if p mod 16 eq 15 then
        return Integers()!((Degree(G)-1)/2);
    else
        return Integers()!(Degree(G)/2);
    end if;
end function;

p := 13;
while p le 10000 do
    p := NextPrime(p);
    num := SspZ8quintic(p);
    if num ne 0 then
        assert p mod 8 eq 7;
        p,num;
    end if;
end while;
quit;
