--
-- See wiki for details on SMART stats: https://en.wikipedia.org/wiki/S.M.A.R.T.
--

DROP TABLE IF EXISTS drive_stats_survival;

CREATE TABLE drive_stats_survival (
    date TEXT NOT NULL,
    serial_number TEXT NOT NULL,
    model TEXT NOT NULL,
    capacity_bytes_MB INTEGER (4) NOT NULL,
	n_fails INTEGER(1) NOT NULL,
	
	is_first_record BOOLEAN NOT NULL,
	is_last_record BOOLEAN NOT NULL,
	min_date TEXT NOT NULL,
    max_date TEXT NOT NULL,
	n_records INTEGER,
    min_hours INTEGER,
    max_hours INTEGER,
	
--	smart_1_raw INTEGER,   -- Read Error Rate
	
	smart_5_raw INTEGER,   -- Reallocated Sectors Count
	
	smart_9_raw INTEGER,   -- Power-On Hours
	
	smart_10_raw INTEGER,  -- Spin Retry Count
	
	smart_12_raw INTEGER,  -- Power Cycle Count
	
	smart_187_raw INTEGER, -- Reported Uncorrectable Errors
	
	smart_188_raw INTEGER, -- Command Timeout
	
	smart_189_raw INTEGER, -- High Fly Writes
	
	smart_196_raw INTEGER, -- Reallocation Event Count
	
	smart_197_raw INTEGER, -- Current Pending Sector Count
	
	smart_198_raw INTEGER, -- (Offline) Uncorrectable Sector Count
	
	smart_201_raw INTEGER, -- Soft Read Error Rate or TA Counter Detected
	
	PRIMARY KEY(date, serial_number)
);	

INSERT INTO drive_stats_survival
	SELECT min(a.date) AS DATE, a.serial_number, a.model, 
	CAST((a.capacity_bytes * 1.0 / 1000000.0) + 0.5 AS INT) AS capacity_bytes_MB,

	b.n_fails,	
	(b.min_date == a.date) AS is_first_record,
	(b.max_date == a.date) AS is_last_record,
	b.min_date,
	b.max_date,
	b.n_records,
	b.min_hours,
	b.max_hours,
		
--	a.smart_1_raw,
	a.smart_5_raw, 
	min(a.smart_9_raw) AS smart_9_raw, 
	a.smart_10_raw,
	a.smart_12_raw, 
	a.smart_187_raw, 
	a.smart_188_raw, 
	a.smart_189_raw,
	a.smart_196_raw,
	a.smart_197_raw, 
	a.smart_198_raw,
	a.smart_201_raw
	
	FROM drive_stats AS 'a'
	INNER JOIN (SELECT 
		c.serial_number,
		sum(c.failure) as n_fails,
		min(c.date) AS min_date,
		max(c.date) AS max_date,
		count(c.date) AS n_records,
		min(c.smart_9_raw) AS min_hours,
		max(c.smart_9_raw) AS max_hours
		FROM drive_stats as 'c'
		GROUP BY c.serial_number) AS 'b'
		ON a.serial_number == b.serial_number
	WHERE 
--		a.serial_number IN ('MJ1311YNG5HWYA', 'W1F0A11P', 'W1F0SCL1', 'STF605MH1MR9DW', '13H2B97AS') AND
		(b.min_date == a.date OR (b.max_date == a.date) OR
			(a.smart_5_raw   IS NOT NULL OR a.smart_5_raw > 0   OR a.smart_5_raw   <> '') OR
			(a.smart_10_raw  IS NOT NULL OR a.smart_10_raw > 0  OR a.smart_10_raw  <> '') OR
--   smart_12_raw is the number of times it has been turned on. Seems to be missing to start with
			(a.smart_187_raw IS NOT NULL OR a.smart_187_raw > 0 OR a.smart_187_raw <> '') OR 
			(a.smart_188_raw IS NOT NULL OR a.smart_188_raw > 0 OR a.smart_188_raw <> '') OR 
			(a.smart_189_raw IS NOT NULL OR a.smart_189_raw > 0 OR a.smart_189_raw <> '') OR
			(a.smart_197_raw IS NOT NULL OR a.smart_197_raw > 0 OR a.smart_197_raw <> '') OR 
			(a.smart_198_raw IS NOT NULL OR a.smart_198_raw > 0 OR a.smart_198_raw <> '') OR
			(a.smart_201_raw IS NOT NULL OR a.smart_201_raw > 0 OR a.smart_201_raw <> ''))
	GROUP BY 
		a.serial_number, 
--		a.smart_1_raw, 
		a.smart_5_raw, 
--		a.smart_9_raw,  -- Dont't want smart 9 - it is the stop time. We take the min value
		a.smart_10_raw,
		a.smart_12_raw, 
		a.smart_187_raw,
		a.smart_188_raw,
		a.smart_189_raw,
		a.smart_198_raw,
		a.smart_201_raw
	ORDER BY a.serial_number, date DESC;
