def stratified_split_by_row(X, y, df, test_size=0.2, random_state=42):
    """Row별로 균형잡힌 분할을 수행하는 함수"""
    X_train, X_test = pd.DataFrame(), pd.DataFrame()
    y_train, y_test = pd.Series(), pd.Series()
    
    for row in ['3d', '4d', '5d']:
        row_mask = df['row'] == row
        X_row = X[row_mask]
        y_row = y[row_mask]
        
        if len(X_row) > 0:
            X_row_train, X_row_test, y_row_train, y_row_test = train_test_split(
                X_row, y_row, 
                test_size=test_size,
                random_state=random_state
            )
            
            X_train = pd.concat([X_train, X_row_train])
            X_test = pd.concat([X_test, X_row_test])
            y_train = pd.concat([y_train, y_row_train])
            y_test = pd.concat([y_test, y_row_test])
    
    return X_train, X_test, y_train, y_test

# main 함수 내에서 데이터 분할 부분을 다음과 같이 수정
X_train, X_test, y_train, y_test = stratified_split_by_row(
    X_imputed, y, df, 
    test_size=args.test_size, 
    random_state=args.random_state
) 