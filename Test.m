function f = Tot(t)
    f = integral(@(x)P(x,t),0,1);
end