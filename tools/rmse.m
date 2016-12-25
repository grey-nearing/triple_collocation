function rmse = rmse(X,Y)

rmse = sqrt((X(:)-Y(:))'*(X(:)-Y(:)));

