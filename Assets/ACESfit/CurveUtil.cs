using System;
using System.Collections.Generic;
using UnityEngine;

#region CurveSegment
public abstract class CurveSegment<T>
{
    public abstract T Eval(T x);
}

public class PowerCurveSegment : CurveSegment<float>
{
    public float m_offsetX;
    public float m_offsetY;
    public float m_scaleX; // always 1 or -1
    public float m_scaleY;
    public float m_lnA;
    public float m_B;

    public float m_ctrPoint;                // Range from 0 to 1
    public float m_ctrPointIntensity;       // Range from -1 to 1

    public PowerCurveSegment()
    {
        Reset();
    }

    public void Reset()
    {
        m_offsetX = 0.0f;
        m_offsetY = 0.0f;
        m_scaleX = 1.0f; // always 1 or -1
        m_lnA = 0.0f;
        m_B = 1.0f;

        m_ctrPoint = 0.5f;
        m_ctrPointIntensity = 0.0f;
    }

    protected float EvalPow(float x)
    {
        float x0 = (x - m_offsetX) * m_scaleX;
        float y0 = 0.0f;

        // log(0) is undefined but our function should evaluate to 0. There are better ways to handle this,
        // but it's doing it the slow way here for clarity.
        if (x0 > 0)
        {
            y0 = Mathf.Exp(m_lnA + m_B * Mathf.Log(x0));
        }

        return y0 * m_scaleY + m_offsetY;
    }

    public override float Eval(float x)
    {
        return EvalPow(x);
    }

    // find a function of the form:
    //   f(x) = e^(lnA + Bln(x))
    // where
    //   f(0)   = 0; not really a constraint
    //   f(x0)  = y0
    //   f'(x0) = m
    public static void SolveAB(ref float lnA, ref float B, float x0, float y0, float m)
    {
        B = (m * x0) / y0;
        lnA = (Mathf.Log(y0)) - (B * Mathf.Log(x0));
    }
}

public class LinearCurveSegment:CurveSegment<float>
{
    public float slope = 1.0f;
    public float intersect = 0.0f;

    public override float Eval(float x)
    {
        return x * slope +intersect;
    }

    public static void AsSlopeIntercept(ref float m, ref float b, float x0, float x1, float y0, float y1)
    {
        float dy = (y1 - y0);
        float dx = (x1 - x0);
        if (dx == 0)
            m = 1.0f;
        else
            m = dy / dx;

        b = y0 - x0 * m;
    }
}

public class ConstantCurveSegment : CurveSegment<float>
{
    public float constValue = 0.0f;
    public override float Eval(float x)
    {
        return constValue;
    }
}

public class CubicBezier : CurveSegment<float>
{
    public Vector2 P0 { set { m_p0 = value; Reconstruct(); } get { return m_p0; } }
    public Vector2 P1 { set { m_p1 = value; Reconstruct(); } get { return m_p1; } }
    public Vector2 P2 { set { m_p2 = value; Reconstruct(); } get { return m_p2; } }
    public Vector2 P3 { set { m_p3 = value; Reconstruct(); } get { return m_p3; } }

    protected Vector2 m_p0;
    protected Vector2 m_p1;
    protected Vector2 m_p2;
    protected Vector2 m_p3;
    protected float m_segment;

    public CubicBezier()
    {
        m_p0 = new Vector2(0, 0);
        m_p1 = new Vector2(1, 1);
        m_p2 = new Vector2(.5f, 0);
        m_p3 = new Vector2(.5f, 1);
        Reconstruct();
    }

    protected virtual void Reconstruct()
    {
        m_segment = CalSegmentRange();
    }

    protected float CalSegmentRange()
    {
        return m_p3.x - m_p0.x;
    }

    public override float Eval(float x)
    {
        float t = EvalItr(x, P0.x, P3.x);
        float _1_t = 1 - t;
        float _1_tMulT = _1_t * t;
        float _1_T_Sqr = _1_t * _1_t;
        return _1_T_Sqr * _1_t * m_p0.y + 3 * _1_tMulT * (_1_t * m_p1.y + t * m_p2.y) + t * t * t * m_p3.y;
    }

    protected float XtoTx(float x){return (x - m_p0.x) / m_segment;}

    protected float EvalItr(float x, float lBound, float rBound, int maxIterationTime = 100, float epsilon = 0.0000001f)
    {
        float itX = x;
        float newX = Lerp(XtoTx(itX)).x;

        float eps = epsilon;
        int maxItTime = maxIterationTime;
        float _lBound = newX > x ? lBound : x;
        float _rBound = newX > x ? x : rBound;
        while (Math.Abs(newX - x)> eps && maxItTime > 0)
        {
            itX = (_lBound + _rBound) / 2.0f;
            newX = Lerp(XtoTx(itX)).x;
            _lBound = newX > x ? _lBound : itX ;
            _rBound = newX > x ? itX : _rBound ;
            maxItTime--;
        }

        return XtoTx(itX);
    }

    public Vector2 Lerp(float t)
    {
        // solve using De Casteljau's Algorithm
        Vector2 p01 = Vector2.Lerp(P0, P1, t);
        Vector2 p12 = Vector2.Lerp(P1, P2, t);
        Vector2 p23 = Vector2.Lerp(P2, P3, t);
        Vector2 p012 = Vector2.Lerp(p01, p12, t);
        Vector2 p123 = Vector2.Lerp(p12, p23, t);
        return Vector2.Lerp(p012, p123, t);
    }
}

public class MonotoneCubicBezier : CurveSegment<float>
{
    protected Vector2 m_angPoint;
    protected float m_ctrStart;
    protected float m_ctrEnd;
    protected CubicBezier m_curve;

    public Vector2 Start
    {
        get { return m_curve.P0; }
        set { m_curve.P0 = value; Set(m_curve.P0, m_curve.P3, m_angPoint, m_ctrStart, m_ctrEnd); }
    }
    public Vector2 End
    {
        get { return m_curve.P3; }
        set { m_curve.P3 = value; Set(m_curve.P0, m_curve.P3, m_angPoint, m_ctrStart, m_ctrEnd); }
    }
    public Vector2 AnglePoint
    {
        get{return m_angPoint;}
        set
        {
            m_angPoint = value;
            Constrain();
            m_curve.P1 = CalculateTangentPoint(m_curve.P0, m_angPoint, m_ctrStart);
            m_curve.P2 = CalculateTangentPoint(m_curve.P0, m_angPoint, m_ctrEnd);
        }
    }
    public float StartHandle
    {
        get { return m_ctrStart; }
        set
        {
            m_ctrStart = value;
            Constrain();
            m_curve.P1 = CalculateTangentPoint(m_curve.P0, m_angPoint, value);
        }
    }
    public float EndHandle
    {
        get { return m_ctrEnd; }
        set
        {
            m_ctrEnd = value;
            Constrain();
            m_curve.P2 = CalculateTangentPoint(m_curve.P3, m_angPoint, value);
        }
    }

    public MonotoneCubicBezier()
    {
        m_curve = new CubicBezier();
        Set(Vector2.zero, Vector2.one, new Vector2(0.5f, 0));
    }

    protected void Constrain()
    {
        float xMax = Mathf.Max(m_curve.P3.x, m_curve.P0.x);
        float xMin = Mathf.Min(m_curve.P3.x, m_curve.P0.x);
        float yMax = Mathf.Max(m_curve.P3.y, m_curve.P0.y);
        float yMin = Mathf.Min(m_curve.P3.y, m_curve.P0.y);
        m_curve.P3.Set(xMax, m_curve.P3.y);
        m_curve.P0.Set(xMin, m_curve.P0.y);
        m_angPoint.x = Mathf.Clamp(m_angPoint.x, xMin, xMax);
        m_angPoint.y = Mathf.Clamp(m_angPoint.y, yMin, yMax);
        m_ctrStart = Mathf.Clamp01(m_ctrStart);
        m_ctrEnd = Mathf.Clamp01(m_ctrEnd);
    }

    protected Vector2 CalculateTangentPoint(Vector2 from , Vector2 to , float t)
    {
        return Vector2.LerpUnclamped(from, to, t);
    }

    public void Set(Vector2 start, Vector2 end, Vector2 m_anglePoint, float startHandle = .5f, float endHandle = .5f)
    {
        m_angPoint = m_anglePoint;
        m_curve.P0 = start;
        m_curve.P3 = end;
        m_ctrStart = startHandle;
        m_ctrEnd = endHandle;
        Constrain();
        m_curve.P1 = CalculateTangentPoint(m_curve.P0, m_anglePoint, m_ctrStart);
        m_curve.P2 = CalculateTangentPoint(m_curve.P3, m_anglePoint, m_ctrEnd);
    }

    public override float Eval(float x)
    {
       return m_curve.Eval(x);
    }
}
    
#endregion

#region ACES curve
public class ACEStonemappingCurve : CurveSegment<float>
{
    public const float original_a = 2.51f;
    public const float original_b = 0.03f;
    public const float original_c = 2.43f;
    public const float original_d = 0.59f;
    public const float original_e = 0.14f;

    public float a;
    public float b;
    public float c;
    public float d;
    public float e;

    /// <summary>
    /// Find if cx^2 +dx + 1 = 0 in fit range.
    /// </summary>
    /// <param name="a"></param>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <param name="d"></param>
    /// <returns></returns>
    public static bool IsParamValid(double a , double b , double c, double d, double maxX ,double minX =0)
    {
        double rootJudge = d * d - 4 * c;
        bool hasRoot = rootJudge >= 0 ? true : false;
        if (hasRoot)
        {
            double x1 = (-d + Math.Sqrt(rootJudge)) / (2 * c);
            double x2 = (-d - Math.Sqrt(rootJudge)) / (2 * c);
            hasRoot = ((x1 > minX && x1 < maxX) || (x2 > minX && x2 < maxX));
        }
        return !hasRoot;
    }

    public bool IsParamValid(double maxX, double minX)
    {
        return IsParamValid(a/e, b/e, c/e, d/e, maxX, minX);
    }

    public ACEStonemappingCurve()
    {
        Reset();
    }

    public void Reset()
    {
        a = original_a;
        b = original_b;
        c = original_c;
        d = original_d;
        e = original_e;
    }

    public override float Eval(float x)
    {
        return ((x * (a * x + b)) / (x * (c * x + d) + e));
    }

    public float CalculateWhitePoint(float epsilon = 1.0f / 255.0f, float maxX = 50.0f, float testStep = 0.01f)
    {
        float x = 0;

        while (this.Eval(x) - 1.0f < epsilon && x < maxX)
        {
            x += testStep;
        }

        return x - testStep;
    }
}
#endregion

#region FilmicCurveOriginal
// ref: http://filmicworlds.com/blog/filmic-tonemapping-with-piecewise-power-curves/
// [Hable17]
// This is how John Hable originally deal with tone mapping curve.
public class FilmicCurveHable : CurveSegment<float>
{
    public Vector2 Toe { get { return m_toe; }set { m_toe = value; Reconstruct(); } }
    public Vector2 Shoulder { get { return m_shoulder; } set { m_shoulder = value; Reconstruct(); } }
    public Vector2 WhitePoint { get { return m_whitePoint; } set { m_whitePoint = value; Reconstruct(); } }
    public float Gamma { get { return m_gamma; } set { m_gamma = value; Reconstruct(); } }

    protected const int m_segmentCount = 4;
    protected Vector2 m_toe, m_shoulder, m_whitePoint;
    protected float m_gamma = 1.0f;

    protected CurveSegment<float>[] m_segments;

    public FilmicCurveHable()
    {
        m_toe = new Vector2(.25f, .25f);
        m_shoulder = new Vector2(.65f, .65f);
        m_whitePoint = new Vector2(1.0f,1.0f);

        m_segments = new CurveSegment<float>[m_segmentCount];
        m_segments[0] = new PowerCurveSegment();
        m_segments[1] = new LinearCurveSegment();
        m_segments[2] = new PowerCurveSegment();
        m_segments[3] = new ConstantCurveSegment();
    }
    public void ConstrainParams()
    {
        m_gamma = Mathf.Max(m_gamma, 0.0f);
    
        // constrain toe into positive
        m_toe.x = Mathf.Max(0.0f, m_toe.x);
        m_toe.y = Mathf.Max(0.0f, m_toe.y);
    
        // constrain shoulder using toe
        m_shoulder.x = Mathf.Max(m_toe.x, m_shoulder.x);
        m_shoulder.y = Mathf.Max(m_toe.y, m_shoulder.y);
    
        // constrain shoulder using toe
        m_whitePoint.x = Mathf.Max(m_whitePoint.x, m_shoulder.x);
        m_whitePoint.y = Mathf.Max(m_whitePoint.y, m_shoulder.y);
    }
    
    public override float Eval(float x)
    {
        int s = x < m_toe.x ? 0 : (x < m_shoulder.x ? 1 : (x < m_whitePoint.x ?  2 : 3));
        return Mathf.Pow( m_segments[s].Eval(x), m_gamma);
    }
    
    protected void Reconstruct()
    {
        Set(m_toe, m_shoulder, m_whitePoint, m_gamma);
    }
    
    public void Set(Vector2 toe, Vector2 shoulder , Vector2 whitePoint, float gamma)
    {
        Set(toe.x, shoulder.x, toe.y, shoulder.y, whitePoint.x, whitePoint.y, gamma);
    }
    
    public void Set(float _x0, float _x1, float _y0, float _y1, float _xw, float _yw = 1.0f, float _gamma = 1.0f)
    {
        PowerCurveSegment toeSeg = m_segments[0] as PowerCurveSegment;
        LinearCurveSegment lSeg = m_segments[1] as LinearCurveSegment;
        PowerCurveSegment shoulderSeg = m_segments[2] as PowerCurveSegment;
        ConstantCurveSegment whitePointSeg = m_segments[3] as ConstantCurveSegment;
    
        //assign params
        {
            this.m_toe.x = _x0;
            this.m_toe.y = _y0;
            this.m_shoulder.x = _x1;
            this.m_shoulder.y = _y1;
            this.m_whitePoint.x = _xw;
            this.m_whitePoint.y = _yw;
            this.m_gamma = _gamma;
            ConstrainParams();
        }

        // linear part
        LinearCurveSegment.AsSlopeIntercept(ref lSeg.slope, ref lSeg.intersect, m_toe.x, m_shoulder.x, m_toe.y, m_shoulder.y);
    
        // toe
        PowerCurveSegment.SolveAB(ref toeSeg.m_lnA, ref toeSeg.m_B, m_toe.x, m_toe.y, lSeg.slope);
        toeSeg.m_offsetX = 0;
        toeSeg.m_offsetY = 0;
        toeSeg.m_scaleX = 1.0f;
        toeSeg.m_scaleY = 1.0f;

        // shoulder
        PowerCurveSegment.SolveAB(ref shoulderSeg.m_lnA, ref shoulderSeg.m_B, m_whitePoint.x - m_shoulder.x, m_whitePoint.y - m_shoulder.y, lSeg.slope);
        shoulderSeg.m_offsetX = m_whitePoint.x;
        shoulderSeg.m_offsetY = m_whitePoint.y;
        shoulderSeg.m_scaleX = -1.0f;
        shoulderSeg.m_scaleY = -1.0f;
    
        // constant
        whitePointSeg.constValue = m_whitePoint.y;
    }
}
#endregion

#region FilmicCurve

public class FilmicCurve : CurveSegment<float>
{
    public Vector2 Toe { get { return m_toe; }set { if (value == m_toe) return; m_toe = value; SetDirty(); } }
    public Vector2 Shoulder { get { return m_shoulder; } }
    public Vector2 ToeAnglePoint { get { return m_toeAngPt; } }
    public Vector2 ShoulderAnglePoint { get { return m_ShoulderAngPt; } }
    public Vector2 WhitePoint { get { return m_whitePoint; } set { if (value == WhitePoint) return; m_whitePoint = value; SetDirty(); } }
    public float ShoulderAngle { get { return m_shoulderAngle; } set { if (value == m_shoulderAngle) return; m_shoulderAngle = value; SetDirty(); } }
    public float ShoulderShoot { get { return m_shoulderShoot; } set { if (value == m_shoulderShoot) return; m_shoulderShoot = value; SetDirty(); } }
    public float ToeStrength { get { return m_toeStr; } set { if (value == m_toeStr) return; m_toeStr = value; SetDirty(); } }
    public float ToeLift { get { return m_toeLift; } set { if (value == m_toeLift) return; m_toeLift = value; SetDirty(); } }
    public float ShoulderStrength { get { return m_shoulderStr; } set { if (value == m_shoulderStr) return; m_shoulderStr = value; SetDirty(); } }
    public float ShoulderLift { get { return m_shoulderLift; } set { if (value == m_shoulderLift) return; m_shoulderLift = value; SetDirty(); } }
    public float Gamma { get { return m_gamma; } set { if (value == m_gamma) return; m_gamma = value; } }

    protected Vector2 m_toe, m_whitePoint, m_blackPoint;
    protected float m_shoulderAngle = 45f;
    protected float m_shoulderShoot = 0.1f;
    protected float m_toeStr = 0.5f;
    protected float m_toeLift = 0.5f;
    protected float m_shoulderStr = 0.5f;
    protected float m_shoulderLift = 0.5f;
    protected float m_gamma = 1.0f;

    protected Vector2 m_toeAngPt, m_ShoulderAngPt, m_shoulder;      // should no be changed from outside.

    protected const int m_segmentCount = 4;
    protected CurveSegment<float>[] m_segments;
    protected bool m_isDataDirty = false;

    public FilmicCurve()
    {
        m_blackPoint = Vector2.zero;
        m_toe = new Vector2(.25f,25f);
        m_shoulder = new Vector2(.5f, .5f);
        m_whitePoint = new Vector2(3.0f, 1.0f);
        m_toeAngPt = Vector2.zero;
        m_ShoulderAngPt = Vector2.zero;

        m_segments = new CurveSegment<float>[4];
        m_segments[0] = new MonotoneCubicBezier();
        m_segments[1] = new LinearCurveSegment();
        m_segments[2] = new MonotoneCubicBezier();
        m_segments[3] = new ConstantCurveSegment();

        Reconstruct();
    }
    protected void SetDirty()
    {
        m_isDataDirty = true;
    }

    protected void Reconstruct()
    {
        m_toe.x = Mathf.Max(m_toe.x, 0.0f);
        m_toe.y = Mathf.Clamp01(m_toe.y);

        // set shoulder:
        m_shoulderAngle = Mathf.Clamp(m_shoulderAngle, 0.001f, 89.99f);
        float rad = Mathf.Deg2Rad * m_shoulderAngle;
        Vector2 shoulderDir = new Vector2(Mathf.Cos(rad), Mathf.Sin(rad));
        float slope = shoulderDir.y / shoulderDir.x;

        // constrain toe point:
        float toeSlope = m_toe.y / m_toe.x;
        if(slope < toeSlope)
        {
            // make the toe slope less equal to linear part slope, this is good for fitting.
            m_toe.y = slope * m_toe.x;
        }

        m_shoulderShoot = Mathf.Clamp(m_shoulderShoot, 0.0f, (1 - m_toe.y) / shoulderDir.y);
        m_shoulder = m_toe + shoulderDir * m_shoulderShoot;

        // set white point:
        m_whitePoint.y = Mathf.Max(m_shoulder.y, m_whitePoint.y);
        float minWhitepointX = m_shoulder.x + (m_whitePoint.y - m_shoulder.y) / shoulderDir.y * shoulderDir.x;
        m_whitePoint.x = Mathf.Max(m_whitePoint.x, minWhitepointX);

        // clamp lift and strength:
        m_toeLift = Mathf.Clamp01(m_toeLift);
        m_toeStr = Mathf.Clamp01(m_toeStr);
        m_shoulderLift = Mathf.Clamp01(m_shoulderLift);
        m_shoulderStr = Mathf.Clamp01(m_shoulderStr);

        // set toe angle point:
        float intersect = m_toe.y - m_toe.x * slope;
        if(intersect < 0)
            m_toeAngPt.Set(-intersect / slope, 0.0f);
        else
            m_toeAngPt.Set(0, intersect);

        // set shoulder angle point:
        m_ShoulderAngPt = new Vector2(minWhitepointX, m_whitePoint.y);

        // set curve segments:
        MonotoneCubicBezier toeCurve =  m_segments[0] as MonotoneCubicBezier;
        toeCurve.Set(m_blackPoint, m_toe, m_toeAngPt, m_toeStr, m_toeLift);

        LinearCurveSegment linearPart = m_segments[1] as LinearCurveSegment;
        linearPart.intersect = intersect;
        linearPart.slope = slope;

        MonotoneCubicBezier shoulderCurve = m_segments[2] as MonotoneCubicBezier;
        shoulderCurve.Set(m_shoulder, m_whitePoint, m_ShoulderAngPt, m_shoulderLift, m_shoulderStr);

        ConstantCurveSegment constPart = m_segments[3] as ConstantCurveSegment;
        constPart.constValue = m_whitePoint.y;

        m_isDataDirty = false;
    }
    public override float Eval(float x)
    {
        if (m_isDataDirty)
            Reconstruct();

        int s = x < m_toe.x ? 0 : (x < m_shoulder.x ? 1 : (x < m_whitePoint.x ? 2 : 3));
        return Mathf.Pow(m_segments[s].Eval(x), m_gamma);
    }
}
#endregion

#region Fitting
public class CurveDataSet
{
    public double[] x;
    public double[] y;
    public double[] dxdy;
    public double[] dx2dy2;

    public int dataCount { get; protected set; }

    public CurveDataSet(int _count)
    {
        dataCount = _count;
        x = new double[dataCount];
        y = new double[dataCount];
        dxdy = new double[dataCount];
        dx2dy2 = new double[dataCount];
    }
}

public class ACESfitting
{
    [Serializable]
    public class Config
    {
        public Vector2 fitRange = new Vector2(0.0f, 5.0f);
        public int initSampleCount = 1000;
        public int estExtraSampleCount = 40;
        public int maxIterationTime = 0;
        public float drevWeight = 1.0f;
        public static Config Default() { return new Config(); }
    }

    protected static double Pow2_d(double x) { return x * x; }
    protected static double Pow3_d(double x) { return x * x * x; }
    protected static double Pow4_d(double x) { double pw2x = x * x; return pw2x * pw2x; }
    protected static double Abs_d(double val) { return (val > 0 ? val : -val); }

    /// <summary>
    /// Calculate Standard Deviation of _curve using data _dataX and _dataY.
    /// </summary>
    /// <param name="_dataX"></param>
    /// <param name="_dataY"></param>
    /// <param name="_curve"></param>
    /// <returns></returns>
    public static double StandardDeviation(List<double> _dataX, List<double> _dataY, CurveSegment<float> _curve)
    {
        if (_dataX == null || _dataY == null || _curve == null || _dataX.Count != _dataY.Count)
            return -1;

        int count = _dataX.Count;
        double sumDevSqr = 0;
        for (int i = 0; i < count; ++i)
        {
            double dev = Abs_d(_curve.Eval((float)_dataX[i]) - _dataY[i]);
            sumDevSqr += dev * dev;
        }
        return Mathf.Sqrt((float)(sumDevSqr / (double)count));
    }

    /// <summary>Computes the solution of a linear equation system.</summary>
    /// <param name="M">
    /// The system of linear equations as an augmented matrix[row, col] where (rows + 1 == cols).
    /// It will contain the solution in "row canonical form" if the function returns "true".
    /// </param>
    /// <returns>Returns whether the matrix has a unique solution or not.</returns>
    protected static bool SolveLinearEquation(double[,] M)
    {
        // input checks
        int rowCount = M.GetUpperBound(0) + 1;
        if (M == null || M.Length != rowCount * (rowCount + 1))
            throw new ArgumentException("The algorithm must be provided with a (n x n+1) matrix.");
        if (rowCount < 1)
            throw new ArgumentException("The matrix must at least have one row.");

        // pivoting
        for (int col = 0; col + 1 < rowCount; col++) if (M[col, col] == 0)
            {
                // find non-zero coefficient
                int swapRow = col + 1;
                for (; swapRow < rowCount; swapRow++) if (M[swapRow, col] != 0) break;

                if (M[swapRow, col] != 0) // found a non-zero coefficient?
                {
                    // yes, then swap it with the above
                    double[] tmp = new double[rowCount + 1];
                    for (int i = 0; i < rowCount + 1; i++)
                    { tmp[i] = M[swapRow, i]; M[swapRow, i] = M[col, i]; M[col, i] = tmp[i]; }
                }
                else return false; // no, then the matrix has no unique solution
            }

        // elimination
        for (int sourceRow = 0; sourceRow + 1 < rowCount; sourceRow++)
        {
            for (int destRow = sourceRow + 1; destRow < rowCount; destRow++)
            {
                double df = M[sourceRow, sourceRow];
                double sf = M[destRow, sourceRow];
                for (int i = 0; i < rowCount + 1; i++)
                {
                    M[destRow, i] = M[destRow, i] * df - M[sourceRow, i] * sf;
                }
            }
        }

        // back-insertion
        for (int row = rowCount - 1; row >= 0; row--)
        {
            double f = M[row, row];
            if (f == 0) return false;

            for (int i = 0; i < rowCount + 1; i++) M[row, i] /= f;
            for (int destRow = 0; destRow < row; destRow++)
            { M[destRow, rowCount] -= M[destRow, row] * M[row, rowCount]; M[destRow, row] = 0; }
        }
        return true;
    }

    // Fit aces curve ((x * (a * x + b)) / (x * (c * x + d) + 1)) with arbitrary data _curveSampleX and _curveSampleY using LSR method.
    // result is stored in result as a,b,c,d.
    protected static bool FitRawACESFunction(List<double> _curveSampleX, List<double> _curveSampleY, out float a, out float b, out float c, out float d)
    {
        // Create linear equation:
        double[,] le = new double[4, 5];
        int sampleCount = _curveSampleX.Count;
        for (int i = 0; i < sampleCount; ++i)
        {
            double xi = _curveSampleX[i], yi = _curveSampleY[i];
            // Please ask razoryang about how the linear equation is constructed.
            // eq1
            le[0, 0] += Pow4_d(xi);         //Sum_pow4x(x, y);
            le[0, 1] += Pow3_d(xi);         // Sum_pow3x(x, y);
            le[0, 2] += -Pow4_d(xi) * yi;   // - Sum_pow4x_y(x, y);
            le[0, 3] += -Pow3_d(xi) * yi;   // - Sum_pow3x_y(_x, _y);
            le[0, 4] += Pow2_d(xi) * yi;    // Sum_pow2x_y(_x, _y);

            // eq2
            le[1, 0] += Pow3_d(xi);         // Sum_pow3x(_x, _y);
            le[1, 1] += Pow2_d(xi);         // Sum_pow2x(_x, _y);
            le[1, 2] += -Pow3_d(xi) * yi;   // -Sum_pow3x_y(_x, _y);
            le[1, 3] += -Pow2_d(xi) * yi;   // - Sum_pow2x_y(_x, _y);
            le[1, 4] += xi * yi;            // Sum_x_y(_x, _y);

            //eq3
            le[2, 0] += Pow4_d(xi) * yi;    // Sum_pow4x_y(_x, _y);
            le[2, 1] += Pow3_d(xi) * yi;    // Sum_pow3x_y(_x, _y);
            le[2, 2] += -Pow4_d(xi) * Pow2_d(yi);     // - Sum_pow4x_pow2y(_x, _y);
            le[2, 3] += -Pow3_d(xi) * Pow2_d(yi);     // - Sum_pow3x_pow2y(_x, _y);
            le[2, 4] += Pow2_d(xi * yi);    // Sum_pow2x_pow2y(_x, _y);

            //eq4
            le[3, 0] += Pow3_d(xi) * yi;    // Sum_pow3x_y(_x, _y);
            le[3, 1] += Pow2_d(xi) * yi;    // Sum_pow2x_y(_x, _y);
            le[3, 2] += -Pow3_d(xi) * Pow2_d(yi);     // - Sum_pow3x_pow2y(_x, _y);
            le[3, 3] += -Pow2_d(xi * yi);       // - Sum_pow2x_pow2y(_x, _y);
            le[3, 4] += xi * Pow2_d(yi);    // Sum_x_pow2y(_x, _y);
        }

        // Solve:
        bool suc = SolveLinearEquation(le);
        a = 0; b =0; c = 0; d = 0;
        if (suc)
        {
            a = (float)le[0, 4];
            b = (float)le[1, 4];
            c = (float)le[2, 4];
            d = (float)le[3, 4];
        }
        return suc;
    }

    // prepare the data for fitting (x, y, dxdy, ddxddy)
    protected static void PrepareCurveData(FilmicCurve _curve, Config _cfg,  ref CurveDataSet _result)
    {
        double _maxX = _cfg.fitRange.y, _minX = _cfg.fitRange.x;
        int _dataCount = _cfg.initSampleCount;
        double range = _maxX - _minX;
        double xStep = range / (double)_dataCount;
        double dyStep = 0.001f;
        for (int i = 0; i < _dataCount; ++i)
        {
            _result.x[i] = _minX + (double)i * xStep;
            _result.y[i] = _curve.Eval((float)_result.x[i]);
            _result.dxdy[i] = (_curve.Eval((float)(_result.x[i] + dyStep)) - _result.y[i]) / dyStep;
        }

        for (int i = 1; i < _dataCount; ++i)
            _result.dx2dy2[i] = (_result.dxdy[i] - _result.dxdy[i - 1]) / xStep;

        _result.dx2dy2[0] = _result.dx2dy2[1];
    }

    protected static void PrepareInitialSampes(CurveSegment<float> _curveFunction, CurveDataSet _curveData, Config _cfg, ref List<double> _resultSampleX, ref List<double> _resultSampleY)
    {
        if (_resultSampleX == null)
            _resultSampleX = new List<double>();
        if (_resultSampleY == null)
            _resultSampleY = new List<double>();

        _resultSampleX.Clear();
        _resultSampleY.Clear();

        int _sampleCnt = _cfg.estExtraSampleCount;
        int xCount = _curveData.dataCount;
        int xCountMinusOne = _curveData.dataCount-1;
        double xStart = _curveData.x[0];
        double xEnd = _curveData.x[xCountMinusOne];
        double xRange = xEnd - xStart;
        double xStep = xRange / xCountMinusOne;

        double d2yWeight = _cfg.drevWeight;
        double exWeight = 0;

        // integral over function
        double sumd2y = 0;
        for (int i = 0; i < xCountMinusOne; ++i)
        {
            sumd2y += 0.5 * Math.Pow((Abs_d(_curveData.dx2dy2[i])* xStep + Abs_d(_curveData.dx2dy2[i+1]) * xStep), d2yWeight) + exWeight;
        }
        int sampleCnt = _sampleCnt;
        
        // calculate per sample area to cover
        double stepArea = sumd2y / (_sampleCnt-1);
        int curSmpIdx = 0;
        int curXptr = 0;
        
        // add sample based on persample covering area:
        while(curSmpIdx < sampleCnt && curXptr < xCount)
        {
            _resultSampleX.Add(_curveData.x[curXptr]);
            _resultSampleY.Add(_curveData.y[curXptr]);

            double xAcc = 0;
            while( xAcc < stepArea && curXptr < (xCount-1) )
            {
                xAcc += 0.5 * Math.Pow((Abs_d(_curveData.dx2dy2[curXptr]) * xStep + Abs_d(_curveData.dx2dy2[curXptr + 1]) * xStep), d2yWeight) + exWeight;
                curXptr++;
            }
            if (curXptr >= xCount - 1)
                curXptr = xCount;
            curSmpIdx++;
        }
    }

    protected static float IterationFit(CurveSegment<float> _originalCurve, CurveDataSet _curveData, List<double> _initialSampleX, List<double> _initialSampleY, Config setting, ref ACEStonemappingCurve _result)
    {
        // initialize:
        List<double> _x = new List<double>(_initialSampleX);
        List<double> _y = new List<double>(_initialSampleY);

        // first Fit
        float a = 0, b = 0, c = 0, d = 0;
        FitRawACESFunction(_x, _y, out a, out b, out c, out d);
        _result.a = a;
        _result.b = b;
        _result.c = c;
        _result.d = d;
        _result.e = 1.0f;

        double finalStdDev = StandardDeviation(new List<double>(_curveData.x), new List<double>(_curveData.y), _result);
        return (float)finalStdDev;
    }

    public static float Fit(FilmicCurve _originalCurve, ref ACEStonemappingCurve _result, Config _cfg, out CurveDataSet _orgData, out List<double> _sampleBufferX , out List<double> _sampleBufferY )
    {
        _sampleBufferX = new List<double>();
        _sampleBufferY = new List<double>();

        _orgData = new CurveDataSet(_cfg.initSampleCount);

        PrepareCurveData(_originalCurve, _cfg, ref _orgData);
        PrepareInitialSampes(_originalCurve, _orgData, _cfg, ref _sampleBufferX, ref _sampleBufferY);
        return IterationFit(_originalCurve, _orgData, _sampleBufferX, _sampleBufferY, _cfg, ref _result); 
    } 

    public static float Fit(FilmicCurve _originalCurve, ref ACEStonemappingCurve _result, Config _cfg)
    {
        CurveDataSet _orgData = null;
        List<double> _sampleBufferX = null;
        List<double> _sampleBufferY = null;

        return Fit(_originalCurve, ref _result, _cfg, out _orgData, out _sampleBufferX, out _sampleBufferY);
    }
}
#endregion

