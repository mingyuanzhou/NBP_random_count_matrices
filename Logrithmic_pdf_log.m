function prob = Logrithmic_pdf_log(u,p)
%prob = 1/(-log(1-p)).*p.^u./u;
prob = -log(-log(1-p))+u.*log(p)-log(u);
