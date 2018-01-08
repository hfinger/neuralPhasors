sigmaHsquare=1;
sigmaVsquare=3;
S_H=4;
S_V=10;
sigmaVHsquare = sigmaVsquare*sigmaHsquare / ( sigmaVsquare + sigmaHsquare );
w_H = (1/sigmaHsquare) / (1/sigmaHsquare+1/sigmaVsquare);
w_V = (1/sigmaVsquare) / (1/sigmaHsquare+1/sigmaVsquare);
S_VH = w_H * S_H + w_V * S_V;

figure(1)
clf;
hold on
plot(gaussian(0:0.1:15,S_H,sqrt(sigmaHsquare)))
plot(gaussian(0:0.1:15,S_V,sqrt(sigmaVsquare)))
plot(gaussian(0:0.1:15,S_VH,sqrt(sigmaVHsquare)))