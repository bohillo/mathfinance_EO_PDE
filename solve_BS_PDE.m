function V = solve_BS_PDE (r0, r1, sigma, term_cond, left_cond, \
			   right_cond, left_bound, right_bound, xmin, xmax, T, Mx, Mt)
	 
	 ## usage: V = solve_BS_PDE (term_cond, left_cond, right_cond,
  ## left_bound, right_bound, T, Mx, Mt)
	 	 
	 V = term_cond';
	 x = linspace(xmin, xmax, Mx + 1)';
	 dt = T / Mt;
	 dx = (xmax - xmin)/Mx;

	 ## Auxillary vectors for numerical scheme
	 is = xmin/dx + (0:Mx)';
	 r = r0 - r1;
	 A = .5*(sigma^2*is.^2-r*is)*dt;
	 B = -(sigma^2*is.^2 + r0)*dt;
	 C = .5*(sigma^2*is.^2 + r*is)*dt;


	 ## Only Crank - Nicholson so far
	 M1 = -spdiags([.5*C(2:Mx), .5*B(2:Mx)-1, .5*A(2:Mx)], [-1 0 1], Mx-1, Mx-1)';
	 M2 = -spdiags([-.5*C(2:Mx), -.5*B(2:Mx)-1, -.5*A(2:Mx)], [-1 0 1], Mx-1, Mx-1)';
	 d = zeros(Mx - 1, 1);

	 for k = Mt:-1:1
	   d(1) = .5*A(1) * (left_cond(k) + left_cond(k+1));
	   d(end) = .5*C(end) * (right_cond(k) + right_cond(k+1));
	   V(2:Mx) = M1\(M2*V(2:Mx) + d);
	 #  V = V .* (x > left_bound(k) - dx / 2 & x < right_bound(k) + dx/2);
	 endfor
					  
	 V = V';

endfunction
