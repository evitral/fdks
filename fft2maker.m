
hx = h0;

hnf = fft2(hx);
emf = hnf.*conj(hnf);
emf2 = emf([129:256 1:128],[129:256 1:128]);

hx = emf2.^0.25;
hy = hx(64:192,64:192);

% surf(emf2.^0.5)