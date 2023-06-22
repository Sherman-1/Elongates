import polars as pl


def create_elongate_database(dataframe : pl.DataFrame):




    filtered = dataframe.filter(
        
        pl.col("type") == "elongate"
        
        
    )