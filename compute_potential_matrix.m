function V = compute_potential_matrix(a,N)
    W = zeros(2*N+1);
    for j=-N:N
        for m=-N:N
            W(N+1+j,N+1+m) = a(2*N+2+j-m);
        end
    end
    W(abs(W)<1e-6) = 0; % Band limit to make V sparse. Remove this for higher accuracy.
    W = sparse(W);
    V = W;

    