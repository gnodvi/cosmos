/* 
 * Файл         : ephaccess.h
 * Исполнитель  : ИПА РАН
 * Заказчик     : ОАО НПК СПП
 * Проект       : КР ФЭЛП, ОКР "Эфемериды"
 * Год          : 2015
 */

#ifndef EPHACCESS_H
#define	EPHACCESS_H

#ifndef EXPORT_SYMBOL
   #ifdef _WIN32
      #define EXPORT_SYMBOL __declspec(dllexport)
   #else
      #define EXPORT_SYMBOL
   #endif
#endif

#ifndef CEXPORT
   #ifndef __cplusplus
      #define CEXPORT EXPORT_SYMBOL
   #else
      #define CEXPORT extern "C" EXPORT_SYMBOL
   #endif
#endif

/* Коды возврата */  
enum
{
  EPH_OK                       =  0,  // успех
  EPH_ERROR_OPENING_FILE       = -1,  // ошибка открытия файла
  EPH_ERROR_READING_FILE       = -2,  // ошибка чтения файла
  EPH_ERROR_BAD_FILE_FORMAT    = -3,  // неподходящий формат файла
  EPH_ERROR_DATE_OUT_OF_RANGE  = -4,  // дата вне допустимого диапазона
  EPH_ERROR_NULL_POINTER       = -5,  // недопустимый нулевой указатель
  EPH_ERROR_BAD_INPUT          = -6,  // недопустимое значение входного параметра
  EPH_ERROR_AMBIGUOUS_CODE     = -7,  // недопустим нулевой идентификатор объекта
  EPH_ERROR_OUT_OF_MEMORY      = -8,  // нехватка оперативной памяти
  EPH_ERROR_UNSUPPORTED_FORMAT = -9,  // неподдерживаемая разновидность формата
  EPH_ERROR_NOT_FOUND          = -10, // в файле не найдена эфемерида требуемого объекта
  EPH_ERROR_INTERNAL_LIMITS    = -11  // превышен внутренний лимит по памяти
};

/* Единицы измерения расстояния и времени */
enum
{
  EPH_AU = 1,
  EPH_KM = 2,
  EPH_SEC = 3,
  EPH_DAY = 4
};

/* Числовые коды небесных тел и других объектов.
 * Нумерация соответствует принятой в формате SPK. */
enum
{
  EPH_SSB        = 0, // Барицентр Солнечной системы
  EPH_MERCURY    = 1, // Меркурий
  EPH_VENUS      = 2, // Венера
  EPH_EMB        = 3, // Барицентр системы Земля-Луна
  EPH_MARS_BC    = 4, // Барицентр системы Марса
  EPH_JUPITER_BC = 5, // Барицентр системы Юпитера
  EPH_SATURN_BC  = 6, // Барицентр системы Сатурна
  EPH_URANUS_BC  = 7, // Барицентр системы Урана
  EPH_NEPTUNE_BC = 8, // Барицентр системы Нептуна
  EPH_PLUTO_BC   = 9, // Барицентр системы Плутона
  EPH_SUN       = 10, // Солнце
  EPH_MOON     = 301, // Луна
  EPH_EARTH    = 399  // Земля
};

/* Числовые коды разниц шкал времени */
enum
{
  EPH_TT_MINUS_TDB = 1000000001 // разность шкал TT - TDB
};

/* Числовые коды лунных систем координат, принятых в различных эфемеридах.
 * Нумерация соответствует соглашениям, принятым производителями
 * эфемерид в формате PCK. */
enum
{
  MOON_PA_DE403   = 31002,
  MOON_PA_DE421   = 31006,
  MOON_PA_DE430   = 32006,
  MOON_PA_INPOP   = 1900301,
  MOON_PA_EPM2011 = 1800301,
  MOON_PA_EPM2015 = 1800302,
  MOON_PA_EPM2017 = 1800303
};

/* Тип структуры, содержащей детали реализации доступа к эфемеридным файлам */
struct tagEphAccess;
typedef struct tagEphAccess EphAccess;

/* ephObjectByName -- определение числового кода объекта по его текстовому
 *                    идентификатору на английском языке.
 * 
 * Входные параметры:
 * name -- C-строка с текстовым идентификатором. NULL не допускается.
 *         Регистр букв не имеет значения. Допустимые идентификаторы:
 *         "SSB", Earth", "Sun", "Moon", "Mars_BC" и т.д., аналогично
 *         числовым константам EPH_SSB и пр.
 * 
 * Возвращаемое значение: числовой код объекта или -1 при отсутствии объекта
 *                        с таким идентификатором.
 */
CEXPORT int ephObjectByName (const char *name);

/* ephCreate -- выделение ресурсов для работы с эфемеридными файлами.
 * 
 * Возвращаемое значение: указатель на структуру EphAccess, содержимое
 *                        которой скрыто для пользователя. При нехватке
 *                        оперативной памяти возвращается NULL.
 * 
 * Примечание: не гарантируется корректная работа при использовании EphAccess *
 *             одновременно в нескольких тредах. В случае необходимости доступа
 *             к эфемеридам параллельно из нескольких тредов следует
 *             вызывать ephCreate для каждого треда в отдельности.
 *             Допускается также вызов ephCreate из одного треда несколько раз
 *             и независимая работа с различными EphAccess * в одном треде.
 */
CEXPORT EphAccess * ephCreate (void);

/* ephLastError -- текст последней ошибки на английском языке
 * 
 * Входные параметры:
 * eph      -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *             NULL не допускается.
 * 
 * Возвращаемое значение: C-строка, содержащая текст ошибки, либо NULL, если
 *                        ошибки при последнем вызове не произошло.
 * Примечание: пользователю не следует освобождать находящуюся по указателю
 *             память.
 */
CEXPORT const char * ephLastError (EphAccess *eph);

/* ephLoadFile -- загрузка эфемеридного файла.
 * 
 * Входные параметры:
 * eph      -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *             NULL не допускается.
 * filename -- путь к файлу, C-строка. NULL не допускается.
 * 
 * Возвращаемое значение: 0 при успешной загрузке, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 */
CEXPORT int ephLoadFile (EphAccess *eph, const char *filename);

/* ephSetDistanceUnits -- установка единиц измерения расстояния.
 * Входные параметры:
 * eph   -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *          NULL не допускается.
 * units -- единицы измерения расстояния. Допускается EPH_KM или EPH_AU,
 *          по умолчанию EPH_KM.
 * Возвращаемое значение: 0 при успешной установке, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 * 
 * Примечание: принято строгое соответствие 1 АЕ = 149597870.7 км,
 * за исключением эфемерид EPM до версии EPM2015, в которых задано
 * собственное значение АЕ.
 */
CEXPORT int ephSetDistanceUnits (EphAccess *eph, int units);

/* ephSetTimeUnits -- установка единиц измерения времени.
 * Входные параметры:
 * eph   -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *          NULL не допускается.
 * units -- единицы измерения времени. Допускается EPH_SEC или EPH_DAY,
 *          по умолчанию EPH_SEC.
 * Возвращаемое значение: 0 при успешной установке, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 */
CEXPORT int ephSetTimeUnits (EphAccess *eph, int units);

/* ephGetLeftmostJD, ephGetRightmostJD -- получение начальной и конечной дат
 *                                        эфемеридного файла.
 * Входные параметры:
 * eph   -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *          NULL не допускается.
 * Возвращаемое значение: значение юлианской даты (-1 в случае, если эфемериды
 * не загружены)
 */
CEXPORT double ephGetLeftmostJD (EphAccess *eph);
CEXPORT double ephGetRightmostJD (EphAccess *eph);

/* ephCalculateRectangular -- вычисление декартовых координат и скоростей
 *                            на заданную дату.
 * Входные параметры:
 * eph       -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *              NULL не допускается.
 * body      -- числовой код искомого объекта
 * reference -- числовой код объекта, относительно которого требуется вычислить
 *              координаты и скорости
 * date0     -- момент времени по шкале TDB в формате юлианской даты (JD),
 *              первое слагаемое
 * date1     -- второе слагаемое момента времени
 * Выходные параметры:
 * xyz       -- массив из трёх компонентов координат (X, Y, Z).
 *              NULL не допускается.
 * vxyz      -- массив из трёх компонентов скоростей (Vx, Vy, Vz).
 *              Допускается NULL, в этом случае скорости не вычисляются.
 * 
 * Примечание:  координаты и скорости искомого тела вычисляются относительно
 *              заданного центра отсчёта, в системе координат, принятой
 *              в используемом эфемеридном файле (как правило, ICRF).
 *              Единицы изменения расстояния и скорости в выходных данных
 *              соответствуют установкам для данного объекта eph, сделанным
 *              в ephSetDistanceUnits и ephSetTimeUnits.
 *              При отсутствии вызова ephSetDistanceUnits расстояние выдаётся
 *              в км. При отсутствии вызова ephSetTimeUnits скорость выдаётся
 *              в км/c.
 * Возвращаемое значение: 0 при успешном вычислении, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 */
CEXPORT int ephCalculateRectangular (EphAccess *eph, int body, int reference,
			     double date0, double date1, double *xyz, double *vxyz);

/* ephCalculateEulerAngles -- вычисление эйлеровых углов и скоростей их
 *                            изменения на заданную дату.
 * Входные параметры:
 * eph   --     указатель на структуру EphAccess, полученный с помощью ephCreate.
 *              NULL не допускается.
 * frame     -- числовой код искомой СК. Допустим 0 в случае, если в загруженном
 *              файле содержится всего одна СК.
 * date0     -- момент времени по шкале TDB в формате юлианской даты (JD),
 *              первое слагаемое
 * date1     -- второе слагаемое момента времени
 * Выходные параметры:
 * angles    -- массив из трёх эйлеровых углов (φ, θ, ψ), определяющих систему
 *              координат относительно родительской (как правило, ICRF).
 *              NULL не допускается. Значения углов выдаются в радианах.
 * dangles   -- массив из трёх скоростей изменения эйлеровых углов.
 *              Допускается NULL, а этом случае скорости изменения углов не
 *              вычисляются.
 *              Значения скоростей изменения выдаются в радианах/день или
 *              радианах/с в зависимости от установок для данного объекта eph,
 *              сделанным в ephSetTimeUnits. При отсутствии вызова
 *              epmSetTimeUnits расстояние скорость выдаётся в радианах/c.
 * Возвращаемое значение: 0 при успешном вычислении, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 */
CEXPORT int ephCalculateEulerAngles (EphAccess *eph, int frame,
            double date0, double date1, double *angles, double *rates);


/* ephCalculateTimeDiff -- вычисление разности шкал времени на заданную дату.
 * Входные параметры:
 * eph       -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *              NULL не допускается.
 * code      -- числовой код искомой разницы шкал. На данный момент допустимо
 *              единственное значение EPH_TT_MINUS_TDB.
 * date0     -- момент времени по шкале TDB в формате юлианской даты (JD),
 *              первое слагаемое
 * date1     -- второе слагаемое момента времени
 * Выходные параметры:
 * diff      -- значение разницы шкал на заданный момент, в единицах измерения,
 *              заданных для данного объекта eph в ephSetTimeUnits.
 *              При отсутствии вызова epmSetTimeUnits разница выдаётся в секундах.
 *              NULL не допускается.
 * Возвращаемое значение: 0 при успешном вычислении, отрицательное значение
 *                        при ошибке (см. константы EPH_ERROR_***)
 */
CEXPORT int ephCalculateTimeDiff (EphAccess *eph, int code, double date0, double date1,
                          double *diff);

/* ephDestroy -- освобождение ресурсов (файловых дескрипторов и памяти).
 * Входные параметры:
 * eph   -- указатель на структуру EphAccess, полученный с помощью ephCreate.
 *          NULL не допускается.
 */
CEXPORT void ephDestroy (EphAccess *eph);

#endif