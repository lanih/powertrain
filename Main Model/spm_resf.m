function Fres = spm_resf(t, y, yp, PTS, iappt)
    Fres = zeros(2*PTS + 2, 1);  % Initialize residual vector
    
    F_const = 96485.;
    % iappt = 27.263836618115; 

    Rp = 2.0e-6;    
    Dsp = 1.0e-14;
    eps_p = 0.385;
    Lp = 80e-6;
    deltar_p = Rp / (PTS-1);
    Ap = (3*(1.0-eps_p))/Rp * Lp;
    iparticle_p = iappt / Ap;


    for i = 1:PTS
        if i == 1
            % HCV at left boundary
            r_face_p = deltar_p / 2;
            term1 = (r_face_p^2 * Dsp) * (y(i+1) - y(i)) / deltar_p;
            term2 = 3 / (r_face_p^3 - 0);
            Fres(i) = term1 * term2 - yp(i);

        elseif i >= 2 && i <= PTS - 1
            % FCV for interior points
            r_p = (i - 1) * deltar_p;
            r_facel_p = r_p - deltar_p / 2;
            r_facer_p = r_p + deltar_p / 2;

            flux_right = r_facer_p^2 * Dsp * (y(i+1) - y(i)) / deltar_p;
            flux_left  = r_facel_p^2 * Dsp * (y(i) - y(i-1)) / deltar_p;

            denom = r_facer_p^3 - r_facel_p^3;
            Fres(i) = (flux_right - flux_left) * (3 / denom) - yp(i);

        else
            % HCV at right boundary
            r_p = (i - 1) * deltar_p;
            r_face_p = r_p - deltar_p / 2;

            term1 = r_p^2 * (iparticle_p / (F_const));
            term2 = r_face_p^2 * Dsp * (y(i) - y(i-1)) / deltar_p;

            denom = r_p^3 - r_face_p^3;
            Fres(i) = (term1 - term2) * (3 / denom) - yp(i);
        end
    end

    Rn = 2.0e-6;    
    Dsn = 2.0e-14;
    eps_n = 0.485;
    Ln = 88e-6;
    deltar_n = Rn / (PTS-1);
    An = (3*(1.0-eps_n))/Rn * Ln;
    iparticle_n = -iappt / An;

    for i = PTS+1:PTS+PTS
        if i == PTS+1
            % HCV at left boundary
            r_face_n = deltar_n / 2;
            term1 = (r_face_n^2 * Dsn) * (y(i+1) - y(i)) / deltar_n;
            term2 = 3 / (r_face_n^3 - 0);
            Fres(i) = term1 * term2 - yp(i);

        elseif i >= PTS+2 && i <= PTS+PTS - 1
            % FCV for interior points
            r_n = (i - PTS - 1) * deltar_n;
            r_facel_n = r_n - deltar_n / 2;
            r_facer_n = r_n + deltar_n / 2;

            flux_right = r_facer_n^2 * Dsn * (y(i+1) - y(i)) / deltar_n;
            flux_left  = r_facel_n^2 * Dsn * (y(i) - y(i-1)) / deltar_n;

            denom = r_facer_n^3 - r_facel_n^3;
            Fres(i) = (flux_right - flux_left) * (3 / denom) - yp(i);

        else
            % HCV at right boundary
            r_n = (i - PTS - 1) * deltar_n;
            r_face_n = r_n - deltar_n / 2;

            term1 = r_n^2 * (iparticle_n / (F_const));
            term2 = r_face_n^2 * Dsn * (y(i) - y(i-1)) / deltar_n;

            denom = r_n^3 - r_face_n^3;
            Fres(i) = (term1 - term2) * (3 / denom) - yp(i);
        end
    end

  POTENTIAL_P = 2*PTS+1;
  POTENTIAL_N = 2*PTS+2;
  RGAS = 8.314;
  T = 298;
    
  % Cathode Anode potential (BV) equation
  KINT_P = 1.04e-11;
  surf_p = y(PTS);
  cmax_p = 51555.;
  j0p = F_const * KINT_P * sqrt(surf_p) * sqrt(cmax_p - surf_p);
  VOLMER_P = (2*RGAS*T / F_const) * asinh(iparticle_p / (2*j0p));
  
  Fres(POTENTIAL_P) = posocp(surf_p / cmax_p) - VOLMER_P - y(POTENTIAL_P);

  % Anode potential (BV) equation
  KINT_N = 2.07e-11;
  surf_n = y(2*PTS);
  cmax_n = 30555.;
  j0n = F_const * KINT_N * sqrt(surf_n) * sqrt(cmax_n - surf_n);
  VOLMER_N = (2*RGAS*T / F_const) * asinh(iparticle_n / (2*j0n));

  Fres(POTENTIAL_N) = negocp(surf_n / cmax_n) - VOLMER_N - y(POTENTIAL_N);

end
