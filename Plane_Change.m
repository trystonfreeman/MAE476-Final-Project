% Calculates requirements for a given planar change
function dv = Plane_Change(mu,r,i_1,i_2,omega_1,omega_2)
    delta = acos(cos(i_1)*cos(i_2)+sin(i_1)*sin(i_2)*cos(omega_1-omega_2));
    v_c = sqrt(mu/r);
    dv = 2*v_c*sin(delta/2);

end