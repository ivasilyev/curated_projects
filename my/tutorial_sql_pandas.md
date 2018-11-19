# A short tutorial for data handling in SQL and `pandas ` 

## 1. Establishing connection

### sql (dBeaver)
```text
Database - Create New Connection - MariaDB - Server Host: genome-euro-mysql.soe.ucsc.edu; Database: hg19; User Name: genomep; Password: password - Next - Next - Finish

SQL Editor - SQL Editor - ';' - Ctrl + Enter
```

### python3
```python
import sqlalchemy
import pandas as pd
import re

engine = sqlalchemy.create_engine('mysql+mysqldb://genomep:password@genome-euro-mysql.soe.ucsc.edu/hg19', pool_recycle=3600)
connection = engine.connect()
```

## 2. Obtaining the whole table
```mysql
SELECT * FROM HInvGeneMrna;
```
```python
df = pd.read_sql('SELECT * FROM HInvGeneMrna;', connection, index_col=None)
```

## 3. Descending sort by matches number:
```mysql
SELECT * FROM HInvGeneMrna ORDER BY matches DESC;
```
```python
df.sort_values(by='matches', ascending=False)
```

## 4. Selecting only rows from chromosome Y:
```mysql
SELECT * FROM HInvGeneMrna WHERE tName = "chrY";
```
```python
df.loc[df['tName'] == 'chrY']
```

## 5. Selecting only rows containing block with size of `163`:
```mysql
SELECT * FROM HInvGeneMrna WHERE blockSizes LIKE "%163%";
```
```python
df.loc[df.loc[:, 'blockSizes'].str.decode(encoding='UTF-8').str.contains('163') == True]
```

## 6. Selecting rows with names ending by `13`:
```mysql
SELECT * FROM HInvGeneMrna WHERE qName REGEXP '.*13$';
```
```python
df.loc[df.loc[:, 'qName'].str.contains('.*13$', regex=True) == True]
```

## 7. Merging tables `HInv` and `HInvGeneMrna` by the column `geneId`/`qName`:
```mysql
SELECT * FROM HInv AS table1
INNER JOIN HInvGeneMrna AS table2 ON table1.geneId = table2.qName;
```
```python
ldf = pd.read_sql('SELECT * FROM HInv;', connection, index_col=None)
rdf = df.copy()
pd.merge(ldf, rdf, left_on='geneId', right_on='qName', how='inner')
```

## 8. Removing the `HIT` substring from the `qName` column:
```mysql
SELECT *, REGEXP_REPLACE(qName, 'HIT', '') FROM HInvGeneMrna;
```
```python
tmp = df.copy()
tmp['qName'] = tmp['qName'].str.replace('HIT', '')
tmp
```

## 9. Selecting only rows with even `mrnaAcc` numeric indexes, merging and nesting queries and aliases:
```mysql
SELECT * FROM
(
	SELECT * FROM
	(
		SELECT *, REGEXP_REPLACE(mrnaAcc, '[^0-9]', '') AS num FROM HInv AS table1 
	) AS table2 
	WHERE num % 2 = 0
) AS table3
INNER JOIN HInvGeneMrna AS table4 ON table3.geneId = table4.qName;
```
```python
pd.merge(ldf.loc[ldf['mrnaAcc'].map(lambda var1: re.findall('\d+', var1)[0]).map(lambda var2: int(var2) % 2 == 0) == True], rdf, left_on='geneId', right_on='qName', how='inner')
```

## 10. Summing the column:
```mysql
SELECT SUM(misMatches) FROM HInvGeneMrna;
```
```python
sum(df.loc[:, 'misMatches'].values.tolist())
```
