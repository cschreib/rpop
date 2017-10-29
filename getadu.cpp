#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 3) {
        print("usage: getadu obs1.fits obs2.fits dark.fits");
        print(" - obs1.fits and obs2.fits must be two repeat observations of the same scene.");
        print("   Ideally this scene will contain a broad range of intensities. Observing");
        print("   parameters must be the same (exposure time, gain, etc.)");
        print(" - dark.fits must be captured with all light sources turned off.");
        print("");
        print("The tool will produce (with bad pixels masked):");
        print(" - sub.fits: 'obs1 - obs2'.");
        print(" - mean.fits: '0.5*(obs1 + obs2)'.");
        print(" - rms.fits: the local pixel RMS in 'sub'.");
        print("");
        print("The measured zero level, ADU and read noise for each chanel will be printed in");
        print("the terminal.");
        return 1;
    }

    fits::input_image a(argv[1]);
    fits::input_image b(argv[2]);
    fits::input_image z(argv[3]);
    fits::output_image r("rms.fits");
    fits::output_image s("sub.fits");
    fits::output_image c("mean.fits");

    for (uint_t i : {1,2,3,4}) {
        a.reach_hdu(i);
        b.reach_hdu(i);
        z.reach_hdu(i);
        r.reach_hdu(i);
        s.reach_hdu(i);
        c.reach_hdu(i);

        // v = gain*adu*n + z
        // where: v = pixel value
        //        n = number of received electrons
        //        z = dark level

        vec2f va, vb, vz;
        a.read(va);
        b.read(vb);
        z.read(vz);

        // Remove dark level with:
        // m = v - z

        double mz = median(vz); {
            vec1u idgz = where(abs(vz - mz) < 20);
            mz = mean(vz[idgz]);
        }

        vec2f m = 0.5*(va + vb) - mz;

        // Poisson + read noise (=rn):
        // sigma(n) = sqrt(n)
        // sigma(m)^2 = (gain*adu)^2*sigma(n)^2 + rn^2
        //            = (gain*adu)^2*n + rn^2 = m*(gain*adu) + rn^2
        // Therefore one can measure 'gain*adu' from the slope of the
        // linear relation between 'm' and 'sigma(m)^2', and the read
        // noise from the intercept.

        // Compute sigma(m) and flag outlier pixels
        vec2f sub = va - vb;

        vec2f rms(sub.dims);
        uint_t hsize = 5;
        vec2f tmp;
        for (uint_t y : range(rms.dims[0]))
        for (uint_t x : range(rms.dims[1])) {
            uint_t x0 = x >= hsize ? x - hsize : 0;
            uint_t x1 = x < rms.dims[1]-hsize ? x + hsize : rms.dims[1]-1;
            uint_t y0 = y >= hsize ? y - hsize : 0;
            uint_t y1 = y < rms.dims[0]-hsize ? y + hsize : rms.dims[0]-1;

            tmp.resize(y1-y0+1, x1-x0+1);
            for (uint_t ty = y0; ty <= y1; ++ty)
            for (uint_t tx = x0; tx <= x1; ++tx) {
                tmp.safe(ty-y0, tx-x0) = sub.safe(ty,tx);
            }

            double v = median(tmp);
            double md = 1.48*median(abs(tmp - v));
            vec1u idg = where(abs(tmp - v) < 5*md);
            rms.safe(y,x) = stddev(tmp[idg]);

            for (uint_t ty = y0; ty <= y1; ++ty)
            for (uint_t tx = x0; tx <= x1; ++tx) {
                tmp.safe(ty-y0, tx-x0) = m.safe(ty,tx);
            }

            v = median(tmp);
            if (abs(m.safe(y,x) - v) >= 4*rms.safe(y,x)) {
                m.safe(y,x) = fnan;
                sub.safe(y,x) = fnan;
            }
        }

        c.write(m);
        s.write(sub);
        r.write(rms);

        // Fit bisector to 'm' vs 'sigma(m)^2'
        vec1u idg = where(is_finite(m));
        auto res1 = linfit(sqr(rms[idg]), 1.0, 1.0, m[idg]);
        auto res2 = linfit(m[idg], 1.0, 1.0, sqr(rms[idg]));

        double s1 = res1.params[1];
        double s2 = 1.0/res2.params[1];
        double o1 = res1.params[0];
        double o2 = res2.params[0];

        double g = (s1*s2 - 1.0 + sqrt((1.0 + sqr(s1))*(1.0 + sqr(s2))))/(s1 + s2);
        double rn = o1 + (o2 + o1/s2)/(1.0 - s1/s2)*(s1 - g);

        // Print data
        print(vec1s{"red", "green1", "green2", "blue"}[i-1]);
        print("dark level: ", mz);
        print("gain x ADU: ", sqrt(g));
        print("read noise: ", sqrt(rn));
    }

    return 0;
}
