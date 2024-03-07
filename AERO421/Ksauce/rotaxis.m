function a = rotaxis(C, phi)
    ax = (C(2,3)-C(3,2))/(2*sin(phi));
    ay = (C(3,1)-C(1,3))/(2*sin(phi));
    az = (C(1,2)-C(2,1))/(2*sin(phi));
    a = [ax ay az]';
end