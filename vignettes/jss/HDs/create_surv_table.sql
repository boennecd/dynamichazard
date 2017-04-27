--
-- See wiki for details on SMART stats: https://en.wikipedia.org/wiki/S.M.A.R.T.
--

DROP TABLE IF EXISTS drive_stats_survival;

CREATE TABLE drive_stats_survival (
    date date NOT NULL,
    serial_number TEXT NOT NULL,
    model TEXT NOT NULL,
    capacity_bytes_GB INTEGER (4) NOT NULL,
	n_fails INTEGER(1) NOT NULL,

	is_first_record BOOLEAN NOT NULL,
	is_last_record BOOLEAN NOT NULL,
	min_date TEXT NOT NULL,
    max_date TEXT NOT NULL,
	n_records INTEGER,
    min_hours INTEGER,
    max_hours INTEGER,

--	smart_1_raw INTEGER,   -- Read Error Rate

  -- Reallocated Sectors Count
  -- From Wiki:
  --    Count of reallocated sectors. The raw value represents a count of the
  --    bad sectors that have been found and remapped.
  -- The count is cummulative so there is no need to aggregate later
	smart_5_raw INTEGER,

	smart_9_raw INTEGER,   -- Power-On Hours

  -- Spin Retry Count
  -- From https://kb.acronis.com/content/9110
  --    S.M.A.R.T. parameter indicates the count of retry of spin start
  --    attempts. This attribute stores a total count of the spin start
  --    attempts to reach the fully operational speed (under the condition that
  --    the first attempt was unsuccessful). Spin attempts are counted for the
  --    entire hard drive's lifetime so far
  -- The count is cummulative so there is no need to aggregate later
	smart_10_raw INTEGER,

	smart_12_raw INTEGER,  -- Power Cycle Count

	smart_187_raw INTEGER, -- Reported Uncorrectable Errors

	smart_188_raw INTEGER, -- Command Timeout

	smart_189_raw INTEGER, -- High Fly Writes:  the count of these errors detected over the lifetime of the drive [i.e. cummulative]


  -- Reallocation Event Count
  -- From wiki:
  --    Count of remap operations. The raw value of this attribute shows the
  --    total count of attempts to transfer data from reallocated sectors to a
  --    spare area. Both successful and unsuccessful attempts are counted.
  -- The count is cummulative so there is no need to aggregate later
	smart_196_raw INTEGER,

	-- Current Pending Sector Count
	-- From Wiki
	--    Count of "unstable" sectors (waiting to be remapped, because of
	--    unrecoverable read errors). If an unstable sector is subsequently read
	--    successfully, the sector is remapped and this value is decreased.
	-- Thus, we always need the present value as the value is decreased and
	-- increased
	smart_197_raw INTEGER,

	smart_198_raw INTEGER, -- (Offline) Uncorrectable Sector Count

	smart_201_raw INTEGER  -- Soft Read Error Rate or TA Counter Detected
);

CREATE INDEX idx_ex101 ON drive_stats(serial_number, date);

INSERT INTO drive_stats_survival
	SELECT
	a.date,
	a.serial_number,
	a.model,
	CAST((a.capacity_bytes * 1.0 / 1000000000.0) + 0.5 AS INT) AS capacity_bytes_MB,

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
	a.smart_9_raw,
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
--		(a.smart_197_raw IS NOT NULL OR a.smart_197_raw > 0 OR a.smart_197_raw <> '') OR
			(a.smart_198_raw IS NOT NULL OR a.smart_198_raw > 0 OR a.smart_198_raw <> '') OR
			(a.smart_201_raw IS NOT NULL OR a.smart_201_raw > 0 OR a.smart_201_raw <> ''))
	ORDER BY a.serial_number, date ASC;

CREATE INDEX idx_ex1 ON drive_stats_survival(serial_number, date);

DROP TABLE IF EXISTS keep_table;

CREATE TEMP TABLE keep_table (
  date date NOT NULL,
  serial_number TEXT NOT NULL);

INSERT INTO keep_table
  SELECT a.date, a.serial_number
  FROM drive_stats_survival AS a
  WHERE NOT EXISTS (
    SELECT *
    FROM (SELECT *
          FROM drive_stats_survival AS c
          WHERE
            c.serial_number == a.serial_number AND
            c.date < a.date
          ORDER BY c.serial_number, c.date DESC
          LIMIT 1) AS b
    WHERE
      b.smart_5_raw   == a.smart_5_raw   AND
      b.smart_10_raw  == a.smart_10_raw  AND
      b.smart_12_raw  == a.smart_12_raw  AND
      b.smart_187_raw == a.smart_187_raw AND
      b.smart_188_raw == a.smart_188_raw AND
      b.smart_189_raw == a.smart_189_raw AND
      b.smart_196_raw == a.smart_196_raw AND
      b.smart_197_raw == a.smart_197_raw AND
      b.smart_198_raw == a.smart_198_raw AND
      b.smart_201_raw == a.smart_201_raw
    );

CREATE INDEX idx_ex2 ON keep_table(serial_number, date);

--SELECT count(*)
--    FROM drive_stats_survival
--	where EXISTS (
--		SELECT *
--		FROM keep_table
--		WHERE
--			drive_stats_survival.serial_number =  keep_table.serial_number AND
--			drive_stats_survival.date =  keep_table.date
--		);

DELETE FROM drive_stats_survival
  WHERE NOT EXISTS (
    SELECT *
    FROM keep_table
    WHERE
      drive_stats_survival.serial_number =  keep_table.serial_number AND
      drive_stats_survival.date =  keep_table.date
    );
