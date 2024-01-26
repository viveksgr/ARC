function [M1_new, M2_new] = ARC_transformMatrix(V, M, b,m)
    % Partition V into V1 (values < 0) and V2 (values >= 0)

    m1 = m(1);
    m2 = m(2);
    m0 = (m1+m2)/2;

    indicesV1 = V < 0;
    indicesV2 = V >= 0;
    V1 = V(indicesV1);
    V2 = V(indicesV2);

    % Corresponding partition of M
    M1 = M(:, indicesV1);
    M2 = M(:, indicesV2);

    M1 = zscore(M1,[],2);
    M2 = zscore(M2,[],2);


    % Bin V1 and V2 and transform M1 and M2
    M1_new = ARC_binAndTransform(M1, V1, b, [m1, m0]);
    M2_new = ARC_binAndTransform(M2, V2, b, [m0, m2]);
end

