function Fres = cathode(t, y, yp)
    % Define or bring in global/constants as needed
    % global PTS deltar_p Dsp iparticle_p n F_const
    % global PTS;
    PTS = 50;

    Rp = 2.0e-6;    
    Dsp = 1.0e-14;
    eps_p = 0.385;
    Lp = 80e-6;
    iappt = 27.263836618115; 
    deltar_p = Rp / (PTS-1);
    Ap = (3*(1.0-eps_p))/Rp * Lp;
    iparticle_p = iappt / Ap;

    F_const = 96485.;

    Fres = zeros(PTS, 1);  % Initialize residual vector

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
end
