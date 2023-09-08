using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Windows;
using System.Collections.Generic;
public class WavesGenerator : MonoBehaviour
{
    public WavesCascade cascade0;
    public WavesCascade cascade1;
    public WavesCascade cascade2;

    // must be a power of 2
    [SerializeField]
    int size = 256;

    [SerializeField]
    WavesSettings wavesSettings;
    [SerializeField]
    bool alwaysRecalculateInitials = false;
    [SerializeField]
    float lengthScale0 = 250;
    [SerializeField]
    float lengthScale1 = 17;
    [SerializeField]
    float lengthScale2 = 5;

    [SerializeField]
    ComputeShader fftShader;
    [SerializeField]
    ComputeShader initialSpectrumShader;
    [SerializeField]
    ComputeShader timeDependentSpectrumShader;
    [SerializeField]
    ComputeShader texturesMergerShader;

    Texture2D gaussianNoise;
    FastFourierTransform fft;
    Texture2D physicsReadback;

    private void Awake()
    {
        Application.targetFrameRate = -1;
        fft = new FastFourierTransform(size, fftShader);
        gaussianNoise = GetNoiseTexture(size);

        cascade0 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);
        cascade1 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);
        cascade2 = new WavesCascade(size, initialSpectrumShader, timeDependentSpectrumShader, texturesMergerShader, fft, gaussianNoise);

        InitialiseCascades();

        physicsReadback = new Texture2D(size, size, TextureFormat.RGBAFloat, false);
    }

    int updateCount = 0;
    void InitialiseCascades()
    {
        updateCount = 0;
        float boundary1 = 2 * Mathf.PI / lengthScale1 * 6f;
        float boundary2 = 2 * Mathf.PI / lengthScale2 * 6f;
        cascade0.CalculateInitials(wavesSettings, lengthScale0, 0.0001f, boundary1);
        cascade1.CalculateInitials(wavesSettings, lengthScale1, boundary1, boundary2);
        cascade2.CalculateInitials(wavesSettings, lengthScale2, boundary2, 9999);

        Shader.SetGlobalFloat("LengthScale0", lengthScale0);
        Shader.SetGlobalFloat("LengthScale1", lengthScale1);
        Shader.SetGlobalFloat("LengthScale2", lengthScale2);
    }

    private void SaveCascadeSequences()
    {
        var sliceCount = 128;
        var current = RenderTexture.active;

        var c0disp = new Texture2D[sliceCount];
        var c1disp = new Texture2D[sliceCount];
        var c2disp = new Texture2D[sliceCount];
        var c0deri = new Texture2D[sliceCount];
        var c1deri = new Texture2D[sliceCount];
        var c2deri = new Texture2D[sliceCount];
        var c0turb = new Texture2D[sliceCount];
        var c1turb = new Texture2D[sliceCount];
        var c2turb = new Texture2D[sliceCount];
        var timeLen = 2.0f;
        float deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade2.CalculateWavesAtTime(time, deltaTime);
            c2deri[i] = CopyRenderTexture(cascade2.Derivatives);
        }

        timeLen = 8.0f;
        deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade0.CalculateWavesAtTime(time, deltaTime);
            cascade1.CalculateWavesAtTime(time, deltaTime);
            cascade2.CalculateWavesAtTime(time, deltaTime);
            c1deri[i] = CopyRenderTexture(cascade1.Derivatives);
        }

        timeLen = 32.0f;
        deltaTime = timeLen / sliceCount;
        for (int i = 0; i < sliceCount; ++i)
        {
            var time = i * deltaTime;
            cascade0.CalculateWavesAtTime(time, deltaTime);
            cascade1.CalculateWavesAtTime(time, deltaTime);
            cascade2.CalculateWavesAtTime(time, deltaTime);

            c0disp[i] = CopyRenderTexture(cascade0.Displacement);
            c0deri[i] = CopyRenderTexture(cascade0.Derivatives);
        }

        //Directory.CreateDirectory("Assets/c0displament");
        BlendCycleTextures(c0disp);
        BlendCycleTextures(c0deri);
        BlendCycleTextures(c1deri);
        BlendCycleTextures(c2deri);

        SaveDerivatives("c0deri", c0deri);
        SaveDerivatives("c1deri", c1deri);
        SaveDerivatives("c2deri", c2deri);
        var count = sliceCount / 2;
        //for (int i = 0; i < count; ++i)
        //{
        //    AssetDatabase.CreateAsset(c0disp[i + count], "Assets/c0displament/t" + i.ToString() + ".asset");
        //}

        var offsets = new List<Color>();
        var scalars = new List<Color>();
        // count must equal 8 * 8
        int columnCount = 8, rowCount = 8;
        var tArray = new Texture2D(size * columnCount, size * rowCount, TextureFormat.RGBAFloat, false);
        var colorData = tArray.GetPixels(0);
        for (int i = 0; i < count; ++i)
        {
            var tex = c0disp[i + count];
            var fData = tex.GetPixels(0);
            var maxColor = new Color(float.MinValue, float.MinValue, float.MinValue, float.MinValue);
            var minColor = new Color(float.MaxValue, float.MaxValue, float.MaxValue, float.MaxValue);
            var colorRange = new Color();

            for (int j = 0; j < fData.Length / 4; ++j)
            {
                var v = fData[j];
                maxColor.r = Mathf.Max(maxColor.r, v.r);
                maxColor.g = Mathf.Max(maxColor.g, v.g);
                maxColor.b = Mathf.Max(maxColor.b, v.b);
                maxColor.a = Mathf.Max(maxColor.a, v.a);

                minColor.r = Mathf.Min(minColor.r, v.r);
                minColor.g = Mathf.Min(minColor.g, v.g);
                minColor.b = Mathf.Min(minColor.b, v.b);
                minColor.a = Mathf.Min(minColor.a, v.a);
            }

            colorRange.r = maxColor.r - minColor.r;
            colorRange.g = maxColor.g - minColor.g;
            colorRange.b = maxColor.b - minColor.b;
            colorRange.a = maxColor.a - minColor.a;

            offsets.Add(minColor);
            scalars.Add(colorRange);

            var column = i % columnCount;
            var row = (rowCount - 1) - i / columnCount;
            for (int j = 0; j < fData.Length; ++j)
            {
                var v = fData[j];
                var r = (v.r - minColor.r) / colorRange.r;
                var g = (v.g - minColor.g) / colorRange.g;
                var b = (v.b - minColor.b) / colorRange.b;
                var a = (v.a - minColor.a) / colorRange.a;

                var localColumn = j % size;
                var localRow = j / size;
                colorData[(row * size + localRow) * (size * columnCount) + (column * size + localColumn)] = v;// new Color(r, g, b, a);
            }
        }
        tArray.SetPixels(colorData);
        tArray.Apply();
        var bytes = tArray.EncodeToEXR();
        System.IO.File.WriteAllBytes("Assets/c0disp.exr", bytes);
        /*var str = new string("{");
        for (int i = 0; i < offsets.Count; ++i)
        {
            var ofs = offsets[i];
            var scalar = scalars[i];
            str += "float4(" + ofs.r.ToString() + ", " + ofs.g.ToString() + ", " + ofs.b.ToString() + ", " + ofs.a.ToString() + ")";
            if (i != offsets.Count - 1)
            {
                str += ", ";
            }
        }
        str += "}";
        Debug.Log(str);*/
        RenderTexture.active = current;
        AssetDatabase.SaveAssets();
        AssetDatabase.Refresh();
    }

    private void SaveDerivatives(string name, Texture2D texture)
    {
        var outputTex = new Texture2D(texture.width, texture.height, TextureFormat.ARGB32, false);
        var colorData = outputTex.GetPixels(0); 
        var fData = texture.GetPixels(0);
        for (int i = 0; i < texture.width; ++i)
        {
            for (int j = 0; j < texture.height; ++j)
            {
                var ofs = j * texture.width + i;
                var v = fData[ofs];
                var r = (v.r * 0.5f + 0.5f);
                var g = (v.g * 0.5f + 0.5f);
                var b = (v.b * 0.5f + 0.5f);
                var a = (v.a * 0.5f + 0.5f);
                colorData[ofs] = new Color(r, g, b, a);
            }
        }
        outputTex.SetPixels(colorData);
        outputTex.Apply();
        var bytes = outputTex.EncodeToPNG();
        System.IO.File.WriteAllBytes("Assets/" + name + ".png", bytes);
    }

    private void SaveDerivatives(string name, Texture2D[] textures)
    {
        var count = textures.Length / 2;
        int columnCount = 8, rowCount = 8;
        var tArray = new Texture2D(size * columnCount, size * rowCount, TextureFormat.RGBAFloat, false);
        var colorData = tArray.GetPixels(0);
        for (int i = 0; i < count; ++i)
        {
            var tex = textures[i + count];
            var fData = tex.GetPixels(0);
            var column = i % columnCount;
            var row = (rowCount - 1) - i / columnCount;
            for (int j = 0; j < fData.Length; ++j)
            {
                var v = fData[j];
                var r = (v.r * 0.5f + 0.5f);
                var g = (v.g * 0.5f + 0.5f);
                var b = (v.b * 0.5f + 0.5f);
                var a = (v.a * 0.5f + 0.5f);

                var localColumn = j % size;
                var localRow = j / size;
                colorData[(row * size + localRow) * (size * columnCount) + (column * size + localColumn)] = v;// new Color(r, g, b, a);
            }
        }
        tArray.SetPixels(colorData);
        tArray.Apply();
        var bytes = tArray.EncodeToEXR();
        System.IO.File.WriteAllBytes("Assets/" + name + ".exr", bytes);
    }

    private void BlendCycleTextures(Texture2D [] texArray)
    {
        var half = texArray.Length / 2;
        float interval = 1.0f / (1 + half);
        for (int i = 0; i < half; ++i)
        {
            var t0Data = texArray[i].GetPixelData<float>(0);
            var factor = interval * i;
            var tex = texArray[i + half];
            var tData = tex.GetPixelData<float>(0);
            float hMin = float.MaxValue, hMax = float.MinValue;
            for (int j = 0; j < tData.Length; ++j)
            {
                tData[j] = Mathf.Lerp(tData[j], t0Data[j], factor);
                hMin = Mathf.Min(hMin, tData[j]);
                hMax = Mathf.Max(hMax, tData[j]);
            }
            tex.SetPixelData<float>(tData, 0);
            tex.Apply();
        }
    }

    private Texture2D CopyRenderTexture(RenderTexture rt)
    {
        RenderTexture.active = rt;
        Texture2D result = new Texture2D(rt.width, rt.height, TextureFormat.RGBAFloat, false);
        result.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, 0);
        result.Apply();
        return result;
    }

    private void SaveRenderTexture(string path, RenderTexture rt)
    {
        RenderTexture.active = rt;
        Texture2D result = new Texture2D(rt.width, rt.height, TextureFormat.RGBAFloat, false);
        result.ReadPixels(new Rect(0, 0, rt.width, rt.height), 0, 0);
        result.Apply();
        AssetDatabase.CreateAsset(result, path);
    }

    private void Update()
    {
        if(UnityEngine.Input.GetKey(KeyCode.P))
        {
            SaveCascadeSequences();
        }

        if (alwaysRecalculateInitials)
        {
            InitialiseCascades();
        }

        cascade0.CalculateWavesAtTime(Time.time, Time.deltaTime);
        //if (updateCount == 0)
        {
            cascade1.CalculateWavesAtTime(Time.time, Time.deltaTime);
            cascade2.CalculateWavesAtTime(Time.time, Time.deltaTime);
        }
        ++updateCount;
        RequestReadbacks();
    }
    Texture2D GetNoiseTexture(int size)
    {
        string filename = "GaussianNoiseTexture" + size.ToString() + "x" + size.ToString();
        Texture2D noise = Resources.Load<Texture2D>("GaussianNoiseTextures/" + filename);
        return noise ? noise : GenerateNoiseTexture(size, true);
    }

    Texture2D GenerateNoiseTexture(int size, bool saveIntoAssetFile)
    {
        Texture2D noise = new Texture2D(size, size, TextureFormat.RGFloat, false, true);
        noise.filterMode = FilterMode.Point;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                noise.SetPixel(i, j, new Vector4(NormalRandom(), NormalRandom()));
            }
        }
        noise.Apply();

#if UNITY_EDITOR
        if (saveIntoAssetFile)
        {
            string filename = "GaussianNoiseTexture" + size.ToString() + "x" + size.ToString();
            string path = "Assets/Resources/GaussianNoiseTextures/";
            AssetDatabase.CreateAsset(noise, path + filename + ".asset");
            Debug.Log("Texture \"" + filename + "\" was created at path \"" + path + "\".");
        }
#endif
        return noise;
    }

    float NormalRandom()
    {
        return Mathf.Cos(2 * Mathf.PI * Random.value) * Mathf.Sqrt(-2 * Mathf.Log(Random.value));
    }

    private void OnDestroy()
    {
        cascade0.Dispose();
        cascade1.Dispose();
        cascade2.Dispose();
    }

    void RequestReadbacks()
    {
        AsyncGPUReadback.Request(cascade0.Displacement, 0, TextureFormat.RGBAFloat, OnCompleteReadback);
    }

    public float GetWaterHeight(Vector3 position)
    {
        Vector3 displacement = GetWaterDisplacement(position);
        displacement = GetWaterDisplacement(position - displacement);
        displacement = GetWaterDisplacement(position - displacement);

        return GetWaterDisplacement(position - displacement).y;
    }

    public Vector3 GetWaterDisplacement(Vector3 position)
    {
        Color c = physicsReadback.GetPixelBilinear(position.x / lengthScale0, position.z / lengthScale0);
        return new Vector3(c.r, c.g, c.b);
    }

    void OnCompleteReadback(AsyncGPUReadbackRequest request) => OnCompleteReadback(request, physicsReadback);

    void OnCompleteReadback(AsyncGPUReadbackRequest request, Texture2D result)
    {
        if (request.hasError)
        {
            Debug.Log("GPU readback error detected.");
            return;
        }
        if (result != null)
        {
            result.LoadRawTextureData(request.GetData<Color>());
            result.Apply();
        }
    }
}
